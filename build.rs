extern crate cc;
/*
extern crate bindgen;


use std::env;
use std::path::PathBuf;
use std::{process::Command, str};
*/
fn main() {

    cc::Build::new()
    .file("src/geodesic.c")
    .compile("libgeodesic.a");
	
	//cc::Build::new()
    //.file("src/TransverseMercatorWrapper.cpp")
    //.compile("libgeotm.a");



/*	
    println!("cargo:rustc-link-lib=llvm");
    println!("cargo:rerun-if-changed=TransverseMercatorWrapper.hpp");

    let llvm_config_out = Command::new("llvm-config")
        .args(&["--cxxflags", "--ldflags", "--system-libs", "--libs", "core"])
        .output()
        .expect("failed to execute llvm-config");

    let llvm_clang_args = llvm_config_out
        .stdout
        .split(|byte| byte.is_ascii_whitespace())
        .map(|arg| str::from_utf8(arg).unwrap());

    let bindings = bindgen::Builder::default()
        .header("TransverseMercatorWrapper.hpp")
        .clang_arg("-x").clang_arg("c++") // c++ flag
        .clang_args(llvm_clang_args)
        .parse_callbacks(Box::new(bindgen::CargoCallbacks))
        .generate()
        .expect("unable to generate bindings");

    let out_path = PathBuf::from(env::var("OUT_DIR").unwrap());
    bindings
        .write_to_file(out_path.join("bindings.rs"))
        .expect("couldn't write bindings!");

*/




}
