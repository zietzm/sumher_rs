fn main() {
    cc::Build::new()
        .file("src/ldak/sumfuns.c")
        .file("src/ldak/decomp.c")
        .compile("ldak");

    println!(r"cargo:rustc-link-lib=static=lapack");
    println!(r"cargo:rustc-link-lib=static=blas");
}
