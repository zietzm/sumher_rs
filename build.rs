fn main() {
    println!("cargo:rerun-if-changed=src/ldak/sumfuns.h");
    println!("cargo:rerun-if-changed=src/ldak/sumfuns.c");
    println!("cargo:rerun-if-changed=src/ldak/decomp.h");
    println!("cargo:rerun-if-changed=src/ldak/decomp.c");
    println!("cargo:rerun-if-changed=src/ldak/fortran_functions.h");

    cc::Build::new()
        .file("src/ldak/sumfuns.c")
        .file("src/ldak/decomp.c")
        .compile("ldak");

    println!(r"cargo:rustc-link-lib=lapack");
    println!(r"cargo:rustc-link-lib=blas");
}
