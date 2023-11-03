use cc;

fn main() {
    cc::Build::new().file("c_src/qfc.c").compile("qfc");
}
