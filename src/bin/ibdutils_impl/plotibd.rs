use axum::extract::Query;
use axum::Json;
use ishare::{
    genome::{self, GenomeInfo},
    gmap,
    indiv::*,
    share::ibd::{ibdseg::IbdSeg, ibdset::*},
};
use serde::{Deserialize, Serialize};

use super::args::*;

pub fn main_plotibd(args: &Commands) {
    if let Commands::PlotIBD {
        genome_info,
        sample_lst1,
        sample_lst2,
        fmt1,
        fmt2,
        ibd1_dir,
        ibd2_dir,
        out,
        sample1,
        sample2,
        port,
    } = args
    {
        use std::sync::Arc;
        let ginfo = Arc::new(genome::GenomeInfo::from_toml_file(&genome_info));
        let gmap = gmap::GeneticMap::from_genome_info(&ginfo);

        let (inds1, inds1_opt) = Individuals::from_txt_file(&sample_lst1);
        let (inds2, inds2_opt) = Individuals::from_txt_file(&sample_lst2);
        let mut ibd1 = IbdSet::new(&gmap, &ginfo, &inds1);
        let mut ibd2 = IbdSet::new(&gmap, &ginfo, &inds2);

        for (((ibd, fmt), dir), inds_opt) in [&mut ibd1, &mut ibd2]
            .iter_mut()
            .zip([fmt1, fmt2])
            .zip([ibd1_dir, ibd2_dir])
            .zip([&inds1_opt, &inds2_opt])
        {
            if fmt.as_str() == "hapibd" {
                ibd.read_hapibd_dir(dir);
                ibd.sort_by_samples();
                ibd.infer_ploidy();
                println!("{:?}", ibd.as_slice().len());
            } else if fmt.as_str() == "tskibd" {
                ibd.read_tskibd_dir(dir);
                ibd.sort_by_haplotypes();
                ibd.infer_ploidy();
            } else if fmt.as_str() == "hmmibd" {
                ibd.read_hmmibd_dir(dir);
                ibd.sort_by_haplotypes();
                ibd.infer_ploidy();
            } else {
                panic!("format {} is not supported.", fmt);
            }
            match inds_opt.as_ref() {
                Some((converter, ind_, PloidConvertDirection::Diploid2Haploid)) => {
                    ibd.covert_to_haploid(ind_, converter);
                }
                Some((converter, ind_, PloidConvertDirection::Haploid2Diploid)) => {
                    ibd.covert_to_het_diploid(ind_, converter);
                }
                None => {}
            }
        }

        assert_eq!(ibd1.get_inds().v(), ibd2.get_inds().v());

        // temporary solution
        let ibdv1 = Arc::new(ibd1.into_vec());
        let ibdv2 = Arc::new(ibd2.into_vec());
        let inds1 = Arc::new(inds1);
        // let inds2 = Arc::new(inds2);

        match port {
            Some(port) => {
                use tokio::runtime::Runtime;
                let rt = Runtime::new().unwrap();
                let url = format!("0.0.0.0:{port}");
                println!("serve at {url}");

                #[derive(Serialize, Deserialize, Debug)]
                struct IdPair {
                    id1: String,
                    id2: String,
                }
                #[derive(Deserialize, Debug)]
                struct SearchQuery {
                    q: String,
                    _page: Option<u32>,
                }

                #[derive(Serialize, Deserialize, Debug)]
                struct Item {
                    id: u32,
                    text: String,
                }
                #[derive(Serialize, Debug)]
                struct SearchResult {
                    items: Vec<Item>,
                    more: bool,
                }

                rt.block_on(async move {
                    use axum::{routing::get, Router};
                    // our router
                    let app = Router::new()
                        .route(
                            "/",
                            get({
                                use axum::response::Html;
                                || async move {
                                    let home =
                                        std::fs::read_to_string("ibdutils_index.html").unwrap();
                                    Html(home)
                                }
                            }),
                        )
                        .route(
                            "/searchid",
                            get({
                                let idv1 = inds1.clone();
                                |query: Query<SearchQuery>| async move {
                                    let q = &query.q;
                                    println!("Query: {:?}", query);
                                    let items = idv1
                                        .as_ref()
                                        .v()
                                        .iter()
                                        .enumerate()
                                        .filter(|(_i, name)| name.contains(q))
                                        .take(5)
                                        .map(|(i, name)| Item {
                                            id: i as u32,
                                            text: name.to_owned(),
                                        })
                                        .collect::<Vec<Item>>();
                                    let search_result = SearchResult {
                                        items,
                                        more: false, // Simulating that there are more pages
                                    };
                                    println!("Results: {:?}", search_result);
                                    Json(search_result)
                                }
                            }),
                        )
                        .route(
                            "/plotibd",
                            get({
                                let ibdv1 = ibdv1.clone();
                                let ibdv2 = ibdv2.clone();
                                let inds1 = inds1.clone();
                                |ids: Query<IdPair>| async move {
                                    let mut id1 = *inds1.m().get(&ids.id1).unwrap_or(&0) as u32;
                                    let mut id2 = *inds1.m().get(&ids.id2).unwrap_or(&1) as u32;
                                    let n = inds1.v().len();
                                    if (id1 == 0) && (id2 == 0) {
                                        use rand::Rng;
                                        let mut rng = rand::thread_rng();
                                        id1 = rng.gen_range(1..n) as u32;
                                        id2 = rng.gen_range(0..id1) as u32;
                                    }
                                    if id1 < id2 {
                                        std::mem::swap(&mut id1, &mut id2);
                                    }
                                    let first = ibdv1[..]
                                        .partition_point(|x| x.individual_pair() < (id1, id2));
                                    let last = ibdv1[..]
                                        .partition_point(|x| x.individual_pair() <= (id1, id2));
                                    let v1 = &ibdv1[first..last];
                                    let first = ibdv2[..]
                                        .partition_point(|x| x.individual_pair() < (id1, id2));
                                    let last = ibdv2[..]
                                        .partition_point(|x| x.individual_pair() <= (id1, id2));
                                    let v2 = &ibdv2[first..last];
                                    let svg_string = plot_svg(v1, v2, ginfo.as_ref()).unwrap();
                                    axum::response::Html(svg_string)
                                }
                            }),
                        );
                    axum::Server::bind(&url.parse().unwrap())
                        .serve(app.into_make_service())
                        .await
                        .unwrap();
                });
            }
            None => {
                let inds = inds1.clone();
                let mut id1 = match sample1.ind_ix1 {
                    Some(id) => {
                        assert!(id < inds.v().len() as u32);
                        id
                    }
                    None => {
                        let s = sample1.ind_name1.as_ref().unwrap();
                        inds.m()[s] as u32
                    }
                };
                let mut id2 = match sample2.ind_ix2 {
                    Some(id) => {
                        assert!(id < inds.v().len() as u32);
                        id
                    }
                    None => {
                        let s = sample2.ind_name2.as_ref().unwrap();
                        inds.m()[s] as u32
                    }
                };
                eprintln!(
                    "Samples\n\tsample1:\t{}\t{}\n\tsample2:\t{}\t{}",
                    inds.v()[id1 as usize],
                    id1,
                    inds.v()[id2 as usize],
                    id2
                );

                if id1 < id2 {
                    std::mem::swap(&mut id1, &mut id2);
                }
                let first = ibdv1[..].partition_point(|x| x.individual_pair() < (id1, id2));
                let last = ibdv1[..].partition_point(|x| x.individual_pair() <= (id1, id2));
                let v1 = &ibdv1[first..last];
                let first = ibdv2[..].partition_point(|x| x.individual_pair() < (id1, id2));
                let last = ibdv2[..].partition_point(|x| x.individual_pair() <= (id1, id2));
                let v2 = &ibdv2[first..last];
                let svg_string = plot_svg(v1, v2, ginfo.as_ref()).unwrap();
                std::fs::write(&out, svg_string)
                    .expect(&format!("cannot write to file: {}", out.to_str().unwrap()));
            }
        }
    }
}
/// plot IBD of a sample pair from both v1 and v2
fn plot_svg(
    v1: &[IbdSeg],
    v2: &[IbdSeg],
    ginfo: &GenomeInfo,
) -> Result<String, Box<dyn std::error::Error>> {
    let mut ret = String::new();

    {
        use plotters::{prelude::*, style::text_anchor::*};
        let root_area = SVGBackend::with_string(&mut ret, (1024, 768)).into_drawing_area();
        root_area.fill(&WHITE)?;

        let nchrom = ginfo.chromnames.len() as f32;
        let chrsz_max = *ginfo.chromsize.iter().max().unwrap() as f32;
        let mut cc = ChartBuilder::on(&root_area)
            .margin(40)
            .set_left_and_bottom_label_area_size(50)
            .set_label_area_size(LabelAreaPosition::Left, 60)
            .caption("IBD shared by a pair of samples", ("sans-serif", 30))
            .build_cartesian_2d(-1.0f32..chrsz_max, -1.0f32..nchrom)?;

        cc.configure_mesh()
            .x_labels(20)
            .y_labels(0)
            .disable_mesh()
            .x_desc("Position")
            .y_desc("Chromosome")
            .axis_desc_style(("sans-serif", 20))
            .x_label_formatter(&|v| format!("{:.0}", v))
            .y_label_formatter(&|v| format!("{:.0}", v))
            .draw()?;

        // draw IBD segments
        let mut points = vec![];
        for (iv, v) in [v1, v2].iter().enumerate() {
            eprintln!("Set{} size: {}", iv, v.len());
            for (iseg, seg) in v.iter().enumerate() {
                let (direction, color, label) = match iv {
                    0 => (1.0f32, &BLUE, "Set1"),
                    1 => (-1.0f32, &RED, "Set2"),
                    _ => panic!(),
                };
                points.clear();
                let (_, m, _, n) = seg.haplotype_pair();
                let (chrid, chrname, s) = ginfo.to_chr_pos(seg.s);
                let e = seg.e - seg.s + s;
                eprintln!("\t{}\t{}\t{}", chrname, s + 1, e + 1);

                let y = chrid as f32 + (m + n + 1) as f32 * 0.03 * direction;
                points.push((s as f32, y));
                points.push((e as f32, y));

                let ls = LineSeries::new(
                    points.iter().map(|x| (x.0, x.1)),
                    color.clone().stroke_width(2),
                );

                let s = cc.draw_series(ls)?;
                let legstyle = |(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], color.clone());
                if iseg == 0 {
                    s.label(label).legend(legstyle);
                }
            }
        }

        let right_center = Pos::new(HPos::Right, VPos::Center);

        // customize y tick labels
        cc.draw_series(PointSeries::of_element(
            (0.0f32..nchrom).step(1.0).values().map(|y| (0.0, y)),
            5,
            ShapeStyle::from(&BLACK).filled(),
            &|coord, size, style| {
                let ts = TextStyle {
                    pos: right_center,
                    font: ("sans-serif", 12).into(),
                    color: RGBAColor(style.color.0, style.color.1, style.color.2, 0.8)
                        .to_backend_color(),
                };
                let chrid = coord.1 as usize;
                let chrname = &ginfo.chromnames[chrid];
                EmptyElement::at(coord)
                    + PathElement::new(vec![(-size, 0), (0, 0)], style)
                    + Text::new(format!("{}", chrname), (-size - 5, 0), ts)
            },
        ))?;

        cc.configure_series_labels().border_style(&BLACK).draw()?;

        // To avoid the IO failure being ignored silently, we manually call the present function
        root_area.present().expect("Unable to write result to file");
    }

    Ok(ret)
}
