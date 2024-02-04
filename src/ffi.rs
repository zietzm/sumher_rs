use crate::io::gwas::GwasSumstats;

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
    gwas_sumstats_1: &GwasSumstats,
    gwas_sumstats_2: &GwasSumstats,
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
            gwas_sumstats_1.sample_size.as_ptr(),
            gwas_sumstats_1.chisq.as_ptr(),
            gwas_sumstats_1.rhos.as_ptr(),
            gwas_sumstats_2.sample_size.as_ptr(),
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_solve_sums_wrapper() {
        let tagging = vec![1.0, 2.0, 3.0];
        let chisqs = vec![1.0, 2.0, 3.0];
        let sample_sizes = vec![1.0, 2.0, 3.0];
        let category_vals = vec![vec![1.0, 2.0, 3.0], vec![4.0, 5.0, 6.0]];
        let category_contribs = vec![vec![1.0, 2.0, 3.0], vec![4.0, 5.0, 6.0]];

        let result = solve_sums_wrapper(
            &tagging,
            &chisqs,
            &sample_sizes,
            &category_vals,
            &category_contribs,
            None,
        )
        .unwrap();

        // Check sizes
        assert_eq!(result.stats.len(), 3 * (2 + 1 + 2));
        assert_eq!(result.likes.len(), 11);
        assert_eq!(result.cohers.len(), 2 * 2);
        assert_eq!(result.influs.len(), 2);

        // Check values
        let desired_stats = [
            0.4791464197764975,
            -0.5588545017977516,
            -0.07970808202125412,
            0.2556046190573969,
            1.3577311773082386,
            2.3675053949121336,
            4.497316861083112,
            2.2773633495918704,
            0.763803940571641,
            5.227708866852336,
            146.47111303693012,
            146.47111303693012,
            0.0,
            87.88266782215808,
            439.41333911079033,
        ];
        for (i, &x) in result.stats.iter().enumerate() {
            assert!((x - desired_stats[i]).abs() < 1e-10);
        }

        let desired_likes = [
            -3.5411094874599045,
            -3.3166019866559435,
            -2.601387310875233,
            0.6760515111840476,
            -0.8075528916021616,
            1.0601008469882118,
            1.0,
            0.0,
            0.0,
            0.0,
            0.0,
        ];
        for (i, &x) in result.likes.iter().enumerate() {
            assert!((x - desired_likes[i]).abs() < 1e-10);
        }

        let desired_cohers = [
            5.605081794938058,
            -10.322278458928105,
            -10.322278458928107,
            20.22585894898246,
        ];
        for (i, &x) in result.cohers.iter().enumerate() {
            assert!((x - desired_cohers[i]).abs() < 1e-10);
        }

        let desired_influs = [4.363636363636364, 1.8545454545454547];
        for (i, &x) in result.influs.iter().enumerate() {
            assert!((x - desired_influs[i]).abs() < 1e-10);
        }
    }
}
