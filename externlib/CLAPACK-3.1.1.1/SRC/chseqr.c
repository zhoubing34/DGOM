#include "f2c.h"
#include "blaswrap.h"

/* Table of constant values */

static complex c_b1 = {0.f,0.f};
static complex c_b2 = {1.f,0.f};
static integer c__1 = 1;
static integer c__2 = 2;
static integer c__49 = 49;

/* Subroutine */ int chseqr_(char *job, char *compz, integer *n, integer *ilo, 
	 integer *ihi, complex *h__, integer *ldh, complex *w, complex *z__, 
	integer *ldz, complex *work, integer *lwork, integer *info)
{
    /* System generated locals */
    address a__1[2];
    integer h_dim1, h_offset, z_dim1, z_offset, i__1, i__2, i__3[2];
    real r__1, r__2, r__3;
    complex q__1;
    char ch__1[2];

    /* Builtin functions */
    /* Subroutine */ int s_cat(char *, char **, integer *, integer *, ftnlen);

    /* Local variables */
    complex hl[2401]	/* was [49][49] */;
    integer kbot, nmin;
    extern logical lsame_(char *, char *);
    extern /* Subroutine */ int ccopy_(integer *, complex *, integer *, 
	    complex *, integer *);
    logical initz;
    complex workl[49];
    logical wantt, wantz;
    extern /* Subroutine */ int claqr0_(logical *, logical *, integer *, 
	    integer *, integer *, complex *, integer *, complex *, integer *, 
	    integer *, complex *, integer *, complex *, integer *, integer *),
	     clahqr_(logical *, logical *, integer *, integer *, integer *, 
	    complex *, integer *, complex *, integer *, integer *, complex *, 
	    integer *, integer *), clacpy_(char *, integer *, integer *, 
	    complex *, integer *, complex *, integer *), claset_(char 
	    *, integer *, integer *, complex *, complex *, complex *, integer 
	    *), xerbla_(char *, integer *);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *);
    logical lquery;


/*  -- LAPACK driver routine (version 3.1) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd.. */
/*     November 2006 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */
/*     Purpose */
/*     ======= */

/*     CHSEQR computes the eigenvalues of a Hessenberg matrix H */
/*     and, optionally, the matrices T and Z from the Schur decomposition */
/*     H = Z T Z**H, where T is an upper triangular matrix (the */
/*     Schur form), and Z is the unitary matrix of Schur vectors. */

/*     Optionally Z may be postmultiplied into an input unitary */
/*     matrix Q so that this routine can give the Schur factorization */
/*     of a matrix A which has been reduced to the Hessenberg form H */
/*     by the unitary matrix Q:  A = Q*H*Q**H = (QZ)*H*(QZ)**H. */

/*     Arguments */
/*     ========= */

/*     JOB   (input) CHARACTER*1 */
/*           = 'E':  compute eigenvalues only; */
/*           = 'S':  compute eigenvalues and the Schur form T. */

/*     COMPZ (input) CHARACTER*1 */
/*           = 'N':  no Schur vectors are computed; */
/*           = 'I':  Z is initialized to the unit matrix and the matrix Z */
/*                   of Schur vectors of H is returned; */
/*           = 'V':  Z must contain an unitary matrix Q on entry, and */
/*                   the product Q*Z is returned. */

/*     N     (input) INTEGER */
/*           The order of the matrix H.  N .GE. 0. */

/*     ILO   (input) INTEGER */
/*     IHI   (input) INTEGER */
/*           It is assumed that H is already upper triangular in rows */
/*           and columns 1:ILO-1 and IHI+1:N. ILO and IHI are normally */
/*           set by a previous call to CGEBAL, and then passed to CGEHRD */
/*           when the matrix output by CGEBAL is reduced to Hessenberg */
/*           form. Otherwise ILO and IHI should be set to 1 and N */
/*           respectively.  If N.GT.0, then 1.LE.ILO.LE.IHI.LE.N. */
/*           If N = 0, then ILO = 1 and IHI = 0. */

/*     H     (input/output) COMPLEX array, dimension (LDH,N) */
/*           On entry, the upper Hessenberg matrix H. */
/*           On exit, if INFO = 0 and JOB = 'S', H contains the upper */
/*           triangular matrix T from the Schur decomposition (the */
/*           Schur form). If INFO = 0 and JOB = 'E', the contents of */
/*           H are unspecified on exit.  (The output value of H when */
/*           INFO.GT.0 is given under the description of INFO below.) */

/*           Unlike earlier versions of CHSEQR, this subroutine may */
/*           explicitly H(i,j) = 0 for i.GT.j and j = 1, 2, ... ILO-1 */
/*           or j = IHI+1, IHI+2, ... N. */

/*     LDH   (input) INTEGER */
/*           The leading dimension of the array H. LDH .GE. max(1,N). */

/*     W        (output) COMPLEX array, dimension (N) */
/*           The computed eigenvalues. If JOB = 'S', the eigenvalues are */
/*           stored in the same order as on the diagonal of the Schur */
/*           form returned in H, with W(i) = H(i,i). */

/*     Z     (input/output) COMPLEX array, dimension (LDZ,N) */
/*           If COMPZ = 'N', Z is not referenced. */
/*           If COMPZ = 'I', on entry Z need not be set and on exit, */
/*           if INFO = 0, Z contains the unitary matrix Z of the Schur */
/*           vectors of H.  If COMPZ = 'V', on entry Z must contain an */
/*           N-by-N matrix Q, which is assumed to be equal to the unit */
/*           matrix except for the submatrix Z(ILO:IHI,ILO:IHI). On exit, */
/*           if INFO = 0, Z contains Q*Z. */
/*           Normally Q is the unitary matrix generated by CUNGHR */
/*           after the call to CGEHRD which formed the Hessenberg matrix */
/*           H. (The output value of Z when INFO.GT.0 is given under */
/*           the description of INFO below.) */

/*     LDZ   (input) INTEGER */
/*           The leading dimension of the array Z.  if COMPZ = 'I' or */
/*           COMPZ = 'V', then LDZ.GE.MAX(1,N).  Otherwize, LDZ.GE.1. */

/*     WORK  (workspace/output) COMPLEX array, dimension (LWORK) */
/*           On exit, if INFO = 0, WORK(1) returns an estimate of */
/*           the optimal value for LWORK. */

/*     LWORK (input) INTEGER */
/*           The dimension of the array WORK.  LWORK .GE. max(1,N) */
/*           is sufficient, but LWORK typically as large as 6*N may */
/*           be required for optimal performance.  A workspace query */
/*           to determine the optimal workspace size is recommended. */

/*           If LWORK = -1, then CHSEQR does a workspace query. */
/*           In this case, CHSEQR checks the input parameters and */
/*           estimates the optimal workspace size for the given */
/*           values of N, ILO and IHI.  The estimate is returned */
/*           in WORK(1).  No error message related to LWORK is */
/*           issued by XERBLA.  Neither H nor Z are accessed. */


/*     INFO  (output) INTEGER */
/*             =  0:  successful exit */
/*           .LT. 0:  if INFO = -i, the i-th argument had an illegal */
/*                    value */
/*           .GT. 0:  if INFO = i, CHSEQR failed to compute all of */
/*                the eigenvalues.  Elements 1:ilo-1 and i+1:n of WR */
/*                and WI contain those eigenvalues which have been */
/*                successfully computed.  (Failures are rare.) */

/*                If INFO .GT. 0 and JOB = 'E', then on exit, the */
/*                remaining unconverged eigenvalues are the eigen- */
/*                values of the upper Hessenberg matrix rows and */
/*                columns ILO through INFO of the final, output */
/*                value of H. */

/*                If INFO .GT. 0 and JOB   = 'S', then on exit */

/*           (*)  (initial value of H)*U  = U*(final value of H) */

/*                where U is a unitary matrix.  The final */
/*                value of  H is upper Hessenberg and triangular in */
/*                rows and columns INFO+1 through IHI. */

/*                If INFO .GT. 0 and COMPZ = 'V', then on exit */

/*                  (final value of Z)  =  (initial value of Z)*U */

/*                where U is the unitary matrix in (*) (regard- */
/*                less of the value of JOB.) */

/*                If INFO .GT. 0 and COMPZ = 'I', then on exit */
/*                      (final value of Z)  = U */
/*                where U is the unitary matrix in (*) (regard- */
/*                less of the value of JOB.) */

/*                If INFO .GT. 0 and COMPZ = 'N', then Z is not */
/*                accessed. */

/*     ================================================================ */
/*             Default values supplied by */
/*             ILAENV(ISPEC,'CHSEQR',JOB(:1)//COMPZ(:1),N,ILO,IHI,LWORK). */
/*             It is suggested that these defaults be adjusted in order */
/*             to attain best performance in each particular */
/*             computational environment. */

/*            ISPEC=1:  The CLAHQR vs CLAQR0 crossover point. */
/*                      Default: 75. (Must be at least 11.) */

/*            ISPEC=2:  Recommended deflation window size. */
/*                      This depends on ILO, IHI and NS.  NS is the */
/*                      number of simultaneous shifts returned */
/*                      by ILAENV(ISPEC=4).  (See ISPEC=4 below.) */
/*                      The default for (IHI-ILO+1).LE.500 is NS. */
/*                      The default for (IHI-ILO+1).GT.500 is 3*NS/2. */

/*            ISPEC=3:  Nibble crossover point. (See ILAENV for */
/*                      details.)  Default: 14% of deflation window */
/*                      size. */

/*            ISPEC=4:  Number of simultaneous shifts, NS, in */
/*                      a multi-shift QR iteration. */

/*                      If IHI-ILO+1 is ... */

/*                      greater than      ...but less    ... the */
/*                      or equal to ...      than        default is */

/*                           1               30          NS -   2(+) */
/*                          30               60          NS -   4(+) */
/*                          60              150          NS =  10(+) */
/*                         150              590          NS =  ** */
/*                         590             3000          NS =  64 */
/*                        3000             6000          NS = 128 */
/*                        6000             infinity      NS = 256 */

/*                  (+)  By default some or all matrices of this order */
/*                       are passed to the implicit double shift routine */
/*                       CLAHQR and NS is ignored.  See ISPEC=1 above */
/*                       and comments in IPARM for details. */

/*                       The asterisks (**) indicate an ad-hoc */
/*                       function of N increasing from 10 to 64. */

/*            ISPEC=5:  Select structured matrix multiply. */
/*                      (See ILAENV for details.) Default: 3. */

/*     ================================================================ */
/*     Based on contributions by */
/*        Karen Braman and Ralph Byers, Department of Mathematics, */
/*        University of Kansas, USA */

/*     ================================================================ */
/*     References: */
/*       K. Braman, R. Byers and R. Mathias, The Multi-Shift QR */
/*       Algorithm Part I: Maintaining Well Focused Shifts, and Level 3 */
/*       Performance, SIAM Journal of Matrix Analysis, volume 23, pages */
/*       929--947, 2002. */

/*       K. Braman, R. Byers and R. Mathias, The Multi-Shift QR */
/*       Algorithm Part II: Aggressive Early Deflation, SIAM Journal */
/*       of Matrix Analysis, volume 23, pages 948--973, 2002. */

/*     ================================================================ */
/*     .. Parameters .. */

/*     ==== Matrices of order NTINY or smaller must be processed by */
/*     .    CLAHQR because of insufficient subdiagonal scratch space. */
/*     .    (This is a hard limit.) ==== */

/*     ==== NL allocates some local workspace to help small matrices */
/*     .    through a rare CLAHQR failure.  NL .GT. NTINY = 11 is */
/*     .    required and NL .LE. NMIN = ILAENV(ISPEC=1,...) is recom- */
/*     .    mended.  (The default value of NMIN is 75.)  Using NL = 49 */
/*     .    allows up to six simultaneous shifts and a 16-by-16 */
/*     .    deflation window.  ==== */

/*     .. */
/*     .. Local Arrays .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     ==== Decode and check the input parameters. ==== */

    /* Parameter adjustments */
    h_dim1 = *ldh;
    h_offset = 1 + h_dim1;
    h__ -= h_offset;
    --w;
    z_dim1 = *ldz;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    --work;

    /* Function Body */
    wantt = lsame_(job, "S");
    initz = lsame_(compz, "I");
    wantz = initz || lsame_(compz, "V");
    r__1 = (real) max(1,*n);
    q__1.r = r__1, q__1.i = 0.f;
    work[1].r = q__1.r, work[1].i = q__1.i;
    lquery = *lwork == -1;

    *info = 0;
    if (! lsame_(job, "E") && ! wantt) {
	*info = -1;
    } else if (! lsame_(compz, "N") && ! wantz) {
	*info = -2;
    } else if (*n < 0) {
	*info = -3;
    } else if (*ilo < 1 || *ilo > max(1,*n)) {
	*info = -4;
    } else if (*ihi < min(*ilo,*n) || *ihi > *n) {
	*info = -5;
    } else if (*ldh < max(1,*n)) {
	*info = -7;
    } else if (*ldz < 1 || wantz && *ldz < max(1,*n)) {
	*info = -10;
    } else if (*lwork < max(1,*n) && ! lquery) {
	*info = -12;
    }

    if (*info != 0) {

/*        ==== Quick return in case of invalid argument. ==== */

	i__1 = -(*info);
	xerbla_("CHSEQR", &i__1);
	return 0;

    } else if (*n == 0) {

/*        ==== Quick return in case N = 0; nothing to do. ==== */

	return 0;

    } else if (lquery) {

/*        ==== Quick return in case of a workspace query ==== */

	claqr0_(&wantt, &wantz, n, ilo, ihi, &h__[h_offset], ldh, &w[1], ilo, 
		ihi, &z__[z_offset], ldz, &work[1], lwork, info);
/*        ==== Ensure reported workspace size is backward-compatible with */
/*        .    previous LAPACK versions. ==== */
/* Computing MAX */
	r__2 = work[1].r, r__3 = (real) max(1,*n);
	r__1 = dmax(r__2,r__3);
	q__1.r = r__1, q__1.i = 0.f;
	work[1].r = q__1.r, work[1].i = q__1.i;
	return 0;

    } else {

/*        ==== copy eigenvalues isolated by CGEBAL ==== */

	if (*ilo > 1) {
	    i__1 = *ilo - 1;
	    i__2 = *ldh + 1;
	    ccopy_(&i__1, &h__[h_offset], &i__2, &w[1], &c__1);
	}
	if (*ihi < *n) {
	    i__1 = *n - *ihi;
	    i__2 = *ldh + 1;
	    ccopy_(&i__1, &h__[*ihi + 1 + (*ihi + 1) * h_dim1], &i__2, &w[*
		    ihi + 1], &c__1);
	}

/*        ==== Initialize Z, if requested ==== */

	if (initz) {
	    claset_("A", n, n, &c_b1, &c_b2, &z__[z_offset], ldz);
	}

/*        ==== Quick return if possible ==== */

	if (*ilo == *ihi) {
	    i__1 = *ilo;
	    i__2 = *ilo + *ilo * h_dim1;
	    w[i__1].r = h__[i__2].r, w[i__1].i = h__[i__2].i;
	    return 0;
	}

/*        ==== CLAHQR/CLAQR0 crossover point ==== */

/* Writing concatenation */
	i__3[0] = 1, a__1[0] = job;
	i__3[1] = 1, a__1[1] = compz;
	s_cat(ch__1, a__1, i__3, &c__2, (ftnlen)2);
	nmin = ilaenv_(&c__1, "CHSEQR", ch__1, n, ilo, ihi, lwork);
	nmin = max(11,nmin);

/*        ==== CLAQR0 for big matrices; CLAHQR for small ones ==== */

	if (*n > nmin) {
	    claqr0_(&wantt, &wantz, n, ilo, ihi, &h__[h_offset], ldh, &w[1], 
		    ilo, ihi, &z__[z_offset], ldz, &work[1], lwork, info);
	} else {

/*           ==== Small matrix ==== */

	    clahqr_(&wantt, &wantz, n, ilo, ihi, &h__[h_offset], ldh, &w[1], 
		    ilo, ihi, &z__[z_offset], ldz, info);

	    if (*info > 0) {

/*              ==== A rare CLAHQR failure!  CLAQR0 sometimes succeeds */
/*              .    when CLAHQR fails. ==== */

		kbot = *info;

		if (*n >= 49) {

/*                 ==== Larger matrices have enough subdiagonal scratch */
/*                 .    space to call CLAQR0 directly. ==== */

		    claqr0_(&wantt, &wantz, n, ilo, &kbot, &h__[h_offset], 
			    ldh, &w[1], ilo, ihi, &z__[z_offset], ldz, &work[
			    1], lwork, info);

		} else {

/*                 ==== Tiny matrices don't have enough subdiagonal */
/*                 .    scratch space to benefit from CLAQR0.  Hence, */
/*                 .    tiny matrices must be copied into a larger */
/*                 .    array before calling CLAQR0. ==== */

		    clacpy_("A", n, n, &h__[h_offset], ldh, hl, &c__49);
		    i__1 = *n + 1 + *n * 49 - 50;
		    hl[i__1].r = 0.f, hl[i__1].i = 0.f;
		    i__1 = 49 - *n;
		    claset_("A", &c__49, &i__1, &c_b1, &c_b1, &hl[(*n + 1) * 
			    49 - 49], &c__49);
		    claqr0_(&wantt, &wantz, &c__49, ilo, &kbot, hl, &c__49, &
			    w[1], ilo, ihi, &z__[z_offset], ldz, workl, &
			    c__49, info);
		    if (wantt || *info != 0) {
			clacpy_("A", n, n, hl, &c__49, &h__[h_offset], ldh);
		    }
		}
	    }
	}

/*        ==== Clear out the trash, if necessary. ==== */

	if ((wantt || *info != 0) && *n > 2) {
	    i__1 = *n - 2;
	    i__2 = *n - 2;
	    claset_("L", &i__1, &i__2, &c_b1, &c_b1, &h__[h_dim1 + 3], ldh);
	}

/*        ==== Ensure reported workspace size is backward-compatible with */
/*        .    previous LAPACK versions. ==== */

/* Computing MAX */
	r__2 = (real) max(1,*n), r__3 = work[1].r;
	r__1 = dmax(r__2,r__3);
	q__1.r = r__1, q__1.i = 0.f;
	work[1].r = q__1.r, work[1].i = q__1.i;
    }

/*     ==== End of CHSEQR ==== */

    return 0;
} /* chseqr_ */
