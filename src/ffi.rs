#[link(name = "ldak")]
extern "C" {
    // void solve_sums(double *stats, double *likes, double *cohers, double *influs,
    //             int num_parts, int gcon, int cept, int num_blocks, int length,
    //             int ncv, int *cvindex, double *cvexps, double *stags,
    //             double **svars, double **ssums, double *snss, double *schis,
    //             int parttype, double tol, int maxiter, int chisol, int sflag,
    //             char *filename)

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
        cvindex: *mut i32,
        cvexps: *mut f64,
        stags: *const f64,
        svars: *const *const f64,
        ssums: *const *const f64,
        snss: *const f64,
        schis: *const f64,
        parttype: i32,
        tol: f64,
        maxiter: i32,
        chisol: i32,
        sflag: i32,
        filename: *const u8,
    );
}

pub struct SolveSumsResult {
    pub stats: Vec<f64>,
    pub likes: Vec<f64>,
    pub cohers: Vec<f64>,
    pub influs: Vec<f64>,
}

pub fn solve_sums_wrapper(
    tagging: &Vec<f64>,
    category_vals: &Vec<Vec<f64>>,
    category_cons: &Vec<Vec<f64>>,
    sample_sizes: &Vec<f64>,
    chisqs: &Vec<f64>,
    n_variants: i32,
    n_categories: i32,
    progress_filename: &str,
    cvindex: Option<&mut Vec<i32>>,
    cvexps: Option<&mut Vec<f64>>,
    gcon: Option<i32>,
    cept: Option<i32>,
    num_blocks: Option<i32>,
    ncv: Option<i32>,
    parttype: Option<i32>,
    tol: Option<f64>,
    maxiter: Option<i32>,
    chisol: Option<i32>,
    sflag: Option<i32>,
) -> SolveSumsResult {
    let mut stats = vec![0.0; 9];
    let mut likes = vec![0.0; 11];
    let mut cohers = vec![0.0; 1];
    let mut influs = vec![0.0; 1];

    let stags_ptr = tagging.as_ptr();
    let svars_ptr = category_vals.as_ptr();
    let ssums_ptr = category_cons.as_ptr();
    let snss_ptr = sample_sizes.as_ptr();
    let schis_ptr = chisqs.as_ptr();
    let cvindex_ptr = match cvindex {
        Some(x) => x.as_mut_ptr(),
        None => std::ptr::null_mut(),
    };
    let cvexps_ptr = match cvexps {
        Some(x) => x.as_mut_ptr(),
        None => std::ptr::null_mut(),
    };

    unsafe {
        solve_sums(
            stats.as_mut_ptr(),
            likes.as_mut_ptr(),
            cohers.as_mut_ptr(),
            influs.as_mut_ptr(),
            n_categories,
            gcon.unwrap_or(0),
            cept.unwrap_or(0),
            num_blocks.unwrap_or(-9999),
            n_variants,
            ncv.unwrap_or(0),
            cvindex_ptr,
            cvexps_ptr,
            stags_ptr,
            svars_ptr as *const *const f64,
            ssums_ptr as *const *const f64,
            snss_ptr,
            schis_ptr,
            parttype.unwrap_or(0),
            tol.unwrap_or(0.001),
            maxiter.unwrap_or(100),
            chisol.unwrap_or(1),
            sflag.unwrap_or(0),
            progress_filename.as_ptr(),
        );
    }

    SolveSumsResult {
        stats,
        likes,
        cohers,
        influs,
    }
}
