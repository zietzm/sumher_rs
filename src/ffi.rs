use crate::io::gwas::AlignedGwasSumstats;

use anyhow::{anyhow, Result};

#[link(name = "ldak")]
extern "C" {
    fn solve_sums(
        stats: *mut f64,
        likes: *mut f64,
        cohers: *mut f64,
        influs: *mut f64,
        num_parts: i32,
        gcon: i32,
        cept: i32,
        num_blocks: i32,
        length: i32,
        ncv: i32,
        cvindex: *const i32,
        cvexps: *const f64,
        stags: *const f64,
        svars: *const *const f64,
        ssums: *const *const f64,
        snss: *const f64,
        schis: *const f64,
        tol: f64,
        maxiter: i32,
        chisol: i32,
        sflag: i32,
    ) -> i32;

    fn solve_cors(
        stats: *mut f64,
        num_parts: i32,
        gcon: i32,
        cept: i32,
        num_blocks: i32,
        length: i32,
        stags: *const f64,
        svars: *const *const f64,
        ssums: *const *const f64,
        snss: *const f64,
        schis: *const f64,
        srhos: *const f64,
        snss2: *const f64,
        schis2: *const f64,
        srhos2: *const f64,
        tol: f64,
        maxiter: i32,
    ) -> i32;
}

fn vec_vec_to_double_ptr_ptr(x: &[Vec<f64>]) -> Vec<*const f64> {
    let mut ptrs = Vec::new();
    for item in x.iter() {
        ptrs.push(item.as_ptr());
    }
    ptrs
}

pub struct SolveSumsOptions {
    pub cvindex: Option<Vec<i32>>,
    pub cvexps: Option<Vec<f64>>,
    pub gcon: Option<i32>,
    pub cept: Option<i32>,
    pub num_blocks: Option<i32>,
    pub ncv: Option<i32>,
    pub tol: Option<f64>,
    pub maxiter: Option<i32>,
    pub chisol: Option<i32>,
    pub sflag: Option<i32>,
}

pub struct SolveCorsOptions {
    pub gcon: Option<i32>,
    pub cept: Option<i32>,
    pub num_blocks: Option<i32>,
    pub parttype: Option<i32>,
    pub tol: Option<f64>,
    pub maxiter: Option<i32>,
}

#[derive(Debug)]
pub struct SolveSumsResult {
    pub stats: Vec<f64>,
    pub likes: Vec<f64>,
    pub cohers: Vec<f64>,
    pub influs: Vec<f64>,
}

#[derive(Debug)]
pub struct SolveCorsResult {
    pub stats: Vec<f64>,
    pub num_parts: i32,
    pub gcon: i32,
    pub cept: i32,
}

pub fn solve_sums_wrapper(
    tagging: &[f64],
    chisqs: &[f64],
    sample_sizes: &[f64],
    category_vals: &[Vec<f64>],
    category_contribs: &[Vec<f64>],
    options: Option<SolveSumsOptions>,
) -> Result<SolveSumsResult> {
    let gcon = options.as_ref().and_then(|x| x.gcon).unwrap_or(0);
    let cept = options.as_ref().and_then(|x| x.cept).unwrap_or(0);
    let num_parts = category_vals.len() as i32;
    let total = num_parts + gcon + cept;

    let mut stats = vec![0.0; 3 * (total + 1 + num_parts) as usize];
    let mut likes = vec![0.0; 11];
    let mut cohers = vec![0.0; num_parts.pow(2) as usize];
    let mut influs = vec![0.0; num_parts as usize];

    let cvindex_ptr = options
        .as_ref()
        .and_then(|x| x.cvindex.as_ref())
        .map_or_else(std::ptr::null, |y| y.as_ptr());

    let cvexps_ptr = options
        .as_ref()
        .and_then(|x| x.cvexps.as_ref())
        .map_or_else(std::ptr::null, |y| y.as_ptr());

    let svars = vec_vec_to_double_ptr_ptr(category_vals);
    let ssums = vec_vec_to_double_ptr_ptr(category_contribs);

    unsafe {
        let code = solve_sums(
            stats.as_mut_ptr(),
            likes.as_mut_ptr(),
            cohers.as_mut_ptr(),
            influs.as_mut_ptr(),
            num_parts,
            gcon,
            cept,
            options.as_ref().and_then(|x| x.num_blocks).unwrap_or(-9999),
            tagging.len() as i32,
            options.as_ref().and_then(|x| x.ncv).unwrap_or(0),
            cvindex_ptr,
            cvexps_ptr,
            tagging.as_ptr(),
            svars.as_ptr(),
            ssums.as_ptr(),
            sample_sizes.as_ptr(),
            chisqs.as_ptr(),
            options.as_ref().and_then(|x| x.tol).unwrap_or(0.001),
            options.as_ref().and_then(|x| x.maxiter).unwrap_or(100),
            options.as_ref().and_then(|x| x.chisol).unwrap_or(1),
            options.as_ref().and_then(|x| x.sflag).unwrap_or(0),
        );
        if code != 0 {
            return Err(anyhow!("solve_sums failed"));
        }
    }

    Ok(SolveSumsResult {
        stats,
        likes,
        cohers,
        influs,
    })
}

pub fn solve_cors_wrapper(
    tagging: &[f64],
    gwas_sumstats_1: &AlignedGwasSumstats,
    gwas_sumstats_2: &AlignedGwasSumstats,
    category_vals: &[Vec<f64>],
    category_contribs: &[Vec<f64>],
    options: Option<SolveCorsOptions>,
) -> Result<SolveCorsResult> {
    let gcon = options.as_ref().and_then(|x| x.gcon).unwrap_or(1);
    let cept = options.as_ref().and_then(|x| x.cept).unwrap_or(0);
    let num_parts = category_vals.len() as i32;
    let total = 2 * (num_parts + gcon + cept) + num_parts + 1;

    let mut stats = vec![0.0; 2 * (total + 4 + num_parts) as usize];

    let svars = vec_vec_to_double_ptr_ptr(category_vals);
    let ssums = vec_vec_to_double_ptr_ptr(category_contribs);

    unsafe {
        let code = solve_cors(
            stats.as_mut_ptr(),
            num_parts,
            gcon,
            cept,
            options.as_ref().and_then(|x| x.num_blocks).unwrap_or(200),
            tagging.len() as i32,
            tagging.as_ptr(),
            svars.as_ptr(),
            ssums.as_ptr(),
            gwas_sumstats_1.sample_sizes.as_ptr(),
            gwas_sumstats_1.chisq.as_ptr(),
            gwas_sumstats_1.rhos.as_ptr(),
            gwas_sumstats_2.sample_sizes.as_ptr(),
            gwas_sumstats_2.chisq.as_ptr(),
            gwas_sumstats_2.rhos.as_ptr(),
            options.as_ref().and_then(|x| x.tol).unwrap_or(0.0001),
            options.as_ref().and_then(|x| x.maxiter).unwrap_or(100),
        );
        if code != 0 {
            return Err(anyhow!("solve_cors failed"));
        }
    }

    Ok(SolveCorsResult {
        stats,
        num_parts,
        gcon,
        cept,
    })
}
