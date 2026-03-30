use error_stack::*;

use super::{IbdUtilsError, Result};

use std::{
    path::{Path, PathBuf},
    sync::Arc,
};

use axum::{
    extract::{Query, State},
    response::{Html, IntoResponse, Response},
    Json,
};

use hyper::StatusCode;
use ishare::{
    genome::{self, GenomeInfo},
    gmap,
    indiv::*,
    share::ibd::{ibdseg::IbdSeg, ibdset::*},
};
use rand::Rng;
use serde::{Deserialize, Serialize};
use tokio::runtime::Runtime;

use super::args::*;

impl IntoResponse for IbdUtilsError {
    fn into_response(self) -> Response {
        let (status, error_message) = match self {
            // Map your errors to status codes here
            IbdUtilsError::ZeroNumChromosome => (StatusCode::BAD_REQUEST, self.to_string()),
            _ => (StatusCode::INTERNAL_SERVER_ERROR, self.to_string()),
        };

        // Return the status and message
        (status, error_message).into_response()
    }
}

pub fn main_plotibd(args: &Commands) -> Result<()> {
    let state = prepare_app_state(args)?;

    match state.port {
        Some(port) => {
            let rt = Runtime::new().change_context(IbdUtilsError::Input)?;
            let url = format!("0.0.0.0:{port}");
            println!("serve at {url}");

            rt.block_on(async move {
                use axum::{routing::get, Router};
                // our router
                let app = Router::new()
                    .route("/", get(handler_html))
                    .route("/searchid", get(handler_searchid))
                    .route("/plotibd", get(handler_plotibd))
                    .with_state(state);

                let listener = tokio::net::TcpListener::bind(&url)
                    .await
                    .change_context(IbdUtilsError::Input)?;
                axum::serve(listener, app)
                    .await
                    .change_context(IbdUtilsError::Input)?;
                Ok::<(), Report<IbdUtilsError>>(())
            })?;
        }
        None => {
            let mut id1 = state.id1;
            let mut id2 = state.id2;
            if id1 < id2 {
                std::mem::swap(&mut id1, &mut id2);
            }
            let first = state.ibd1[..].partition_point(|x| x.individual_pair() < (id1, id2));
            let last = state.ibd1[..].partition_point(|x| x.individual_pair() <= (id1, id2));
            let v1 = &state.ibd1[first..last];
            let first = state.ibd2[..].partition_point(|x| x.individual_pair() < (id1, id2));
            let last = state.ibd2[..].partition_point(|x| x.individual_pair() <= (id1, id2));
            let v2 = &state.ibd2[first..last];
            let svg_string =
                plot_svg(v1, v2, state.ginfo.as_ref()).change_context(IbdUtilsError::Input)?;
            std::fs::write(state.out.clone(), svg_string).change_context(IbdUtilsError::Input)?;
        }
    }

    Ok(())
}

fn prepare_app_state(args: &Commands) -> Result<AppState> {
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
        let ginfo = Arc::new(
            genome::GenomeInfo::from_toml_file(genome_info).change_context(IbdUtilsError::Input)?,
        );
        let gmap = Arc::new(
            gmap::GeneticMap::from_genome_info(&ginfo).change_context(IbdUtilsError::Input)?,
        );

        let (inds1, inds1_opt) =
            Individuals::from_txt_file(sample_lst1).change_context(IbdUtilsError::Input)?;
        let inds1 = Arc::new(inds1);
        let (inds2, inds2_opt) =
            Individuals::from_txt_file(sample_lst2).change_context(IbdUtilsError::Input)?;
        let inds2 = Arc::new(inds2);
        let mut ibd1 = IbdSet::new(gmap.clone(), ginfo.clone(), inds1.clone());
        let mut ibd2 = IbdSet::new(gmap, ginfo.clone(), inds2);

        for (((ibd, fmt), dir), inds_opt) in [&mut ibd1, &mut ibd2]
            .iter_mut()
            .zip([fmt1, fmt2])
            .zip([ibd1_dir, ibd2_dir])
            .zip([inds1_opt, inds2_opt].into_iter())
        {
            if fmt.as_str() == "hapibd" {
                ibd.read_hapibd_dir(dir)
                    .change_context(IbdUtilsError::Input)?;
                ibd.sort_by_samples();
                ibd.infer_ploidy();
            } else if fmt.as_str() == "tskibd" {
                ibd.read_tskibd_dir(dir)
                    .change_context(IbdUtilsError::Input)?;
                ibd.sort_by_haplotypes();
                ibd.infer_ploidy();
            } else if fmt.as_str() == "hmmibd" {
                ibd.read_hmmibd_dir(dir)
                    .change_context(IbdUtilsError::Input)?;
                ibd.sort_by_haplotypes();
                ibd.infer_ploidy();
            } else {
                panic!("format {fmt} is not supported.");
            }
            match inds_opt {
                Some((converter, ind_, PloidConvertDirection::Diploid2Haploid)) => {
                    ibd.covert_to_haploid(Arc::new(ind_), &converter);
                }
                Some((converter, ind_, PloidConvertDirection::Haploid2Diploid)) => {
                    ibd.covert_to_het_diploid(Arc::new(ind_), &converter)
                        .change_context(IbdUtilsError::Library)?;
                }
                None => {}
            }
        }

        ensure!(
            ibd1.get_inds().v() == ibd2.get_inds().v(),
            IbdUtilsError::Input
                .into_report()
                .attach("Ibd1IndNotEqualId2Indvi")
        );
        let inds = inds1.clone();
        let id1 = match sample1.ind_ix1 {
            Some(id) => {
                ensure!(
                    id < inds.v().len() as u32,
                    IbdUtilsError::Input.into_report().attach("InvalidSampleId")
                );
                id
            }
            None => {
                let s = sample1
                    .ind_name1
                    .as_ref()
                    .ok_or(IbdUtilsError::Input)
                    .attach("MissBothSampleIdAndNameSnafu")?;
                *inds
                    .m()
                    .get(s)
                    .ok_or(IbdUtilsError::Input)
                    .attach("InvalidSampleName")? as u32
            }
        };
        let id2 = match sample2.ind_ix2 {
            Some(id) => {
                ensure!(
                    id < inds.v().len() as u32,
                    IbdUtilsError::Input.into_report().attach("InvalidSampleId")
                );
                id
            }
            None => {
                let s = sample2
                    .ind_name2
                    .as_ref()
                    .ok_or(IbdUtilsError::Input)
                    .attach("InvalidSampleIdandName")?;
                *inds
                    .m()
                    .get(s)
                    .ok_or(IbdUtilsError::Input)
                    .attach("InvalidSampleId")? as u32
            }
        };
        eprintln!(
            "Samples\n\tsample1:\t{}\t{}\n\tsample2:\t{}\t{}",
            inds.v()[id1 as usize],
            id1,
            inds.v()[id2 as usize],
            id2
        );

        // temporary solution
        Ok(AppState {
            ibd1: Arc::from(ibd1.into_vec().into_boxed_slice()),
            ibd2: Arc::from(ibd2.into_vec().into_boxed_slice()),
            inds: inds1.clone(),
            port: *port,
            id1,
            id2,
            out: Arc::from(out.to_path_buf().into_boxed_path()),
            ginfo,
        })
    } else {
        // should not happening
        Ok(AppState {
            ibd1: Arc::from(vec![].into_boxed_slice()),
            ibd2: Arc::from(vec![].into_boxed_slice()),
            inds: Arc::new(Individuals::from_str_iter(vec![].into_iter())),
            port: None,
            id1: 0,
            id2: 0,
            out: Arc::from(PathBuf::new().into_boxed_path()),
            ginfo: Arc::new(GenomeInfo::new()),
        })
    }
}

/// plot IBD of a sample pair from both v1 and v2
fn plot_svg(v1: &[IbdSeg], v2: &[IbdSeg], ginfo: &GenomeInfo) -> Result<String> {
    let mut ret = String::new();

    {
        use plotters::{prelude::*, style::text_anchor::*};
        let root_area = SVGBackend::with_string(&mut ret, (1024, 768)).into_drawing_area();
        root_area
            .fill(&WHITE)
            .change_context(IbdUtilsError::Input)
            .attach("error in fill root area")?;

        let nchrom = ginfo.chromnames.len() as f32;
        let chrsz_max = *ginfo
            .chromsize
            .iter()
            .max()
            .ok_or(IbdUtilsError::Library)
            .attach("zero number of chromosome")? as f32;
        let mut cc = ChartBuilder::on(&root_area)
            .margin(40)
            .set_left_and_bottom_label_area_size(50)
            .set_label_area_size(LabelAreaPosition::Left, 60)
            .caption("IBD shared by a pair of samples", ("sans-serif", 30))
            .build_cartesian_2d(-1.0f32..chrsz_max, -1.0f32..nchrom)
            .change_context(IbdUtilsError::Library)
            .attach("error building cartesian 2d")?;

        cc.configure_mesh()
            .x_labels(20)
            .y_labels(0)
            .disable_mesh()
            .x_desc("Position")
            .y_desc("Chromosome")
            .axis_desc_style(("sans-serif", 20))
            .x_label_formatter(&|v| format!("{v:.0}"))
            .y_label_formatter(&|v| format!("{v:.0}"))
            .draw()
            .change_context(IbdUtilsError::Library)
            .attach("draw error")?;

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
                let (id1, m, id2, n) = seg.haplotype_pair();
                let (chrid, chrname, s) = ginfo.to_chr_pos(seg.s);
                let e = seg.e - seg.s + s;
                eprintln!("\t{}\t{}\t{}\t{}\t{}", id1, id2, chrname, s + 1, e + 1);

                let y = chrid as f32 + (m + n + 1) as f32 * 0.03 * direction;
                points.push((s as f32, y));
                points.push((e as f32, y));

                let ls = LineSeries::new(
                    points.iter().map(|x| (x.0, x.1)),
                    color.clone().stroke_width(2),
                );

                let s = cc
                    .draw_series(ls)
                    .change_context(IbdUtilsError::Library)
                    .attach("error in draw_series")?;
                let legstyle = |(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], *color);
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
                    + Text::new(chrname.to_string(), (-size - 5, 0), ts)
            },
        ))
        .change_context(IbdUtilsError::Library)
        .attach("error in draw_series")?;

        cc.configure_series_labels()
            .border_style(BLACK)
            .draw()
            .change_context(IbdUtilsError::Input)
            .attach("error in draw")?;

        // To avoid the IO failure being ignored silently, we manually call the present function
        root_area
            .present()
            .change_context(IbdUtilsError::Input)
            .attach("Unable to write result to file")?;
    }

    Ok(ret)
}

async fn handler_html() -> std::result::Result<Html<String>, String> {
    let home = std::fs::read_to_string("ibdutils_index.html")
        .change_context(IbdUtilsError::Input)
        .map_err(|e| format!("{e:?}"))?;
    Ok(Html(home))
}

async fn handler_searchid(query: Query<SearchQuery>, state: State<AppState>) -> Json<SearchResult> {
    let q = &query.q;
    println!("Query: {query:?}");
    let items = state
        .inds
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
    println!("Results: {search_result:?}");
    Json(search_result)
}
async fn handler_plotibd(
    ids: Query<IdPair>,
    state: State<AppState>,
) -> std::result::Result<axum::response::Html<String>, String> {
    let inds = &state.inds;
    let mut id1: u32 = ids
        .id1
        .parse()
        .change_context(IbdUtilsError::Input)
        .map_err(|e| format!("{e:?}"))?;
    let mut id2: u32 = ids
        .id2
        .parse()
        .change_context(IbdUtilsError::Input)
        .map_err(|e| format!("{e:?}"))?;
    let n = inds.v().len();
    if (id1 == 0) && (id2 == 0) {
        let mut rng = rand::rng();
        id1 = rng.random_range(1..n) as u32;
        id2 = rng.random_range(0..id1) as u32;
    }
    if id1 < id2 {
        std::mem::swap(&mut id1, &mut id2);
    }
    let first = state.ibd1[..].partition_point(|x| x.individual_pair() < (id1, id2));
    let last = state.ibd1[..].partition_point(|x| x.individual_pair() <= (id1, id2));
    let v1 = &state.ibd1[first..last];
    let first = state.ibd2[..].partition_point(|x| x.individual_pair() < (id1, id2));
    let last = state.ibd2[..].partition_point(|x| x.individual_pair() <= (id1, id2));
    let v2 = &state.ibd2[first..last];
    let svg_string = plot_svg(v1, v2, state.ginfo.as_ref())
        .change_context(IbdUtilsError::Input)
        .map_err(|e| format!("{e:?}"))?;
    Ok(axum::response::Html(svg_string))
}

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

#[derive(Clone, Debug)]
struct AppState {
    ibd1: Arc<[IbdSeg]>,
    ibd2: Arc<[IbdSeg]>,
    inds: Arc<Individuals>,
    port: Option<u16>,
    id1: u32,
    id2: u32,
    out: Arc<Path>,
    ginfo: Arc<GenomeInfo>,
}
