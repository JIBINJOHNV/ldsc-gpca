root <- commandArgs(trailingOnly = TRUE)[1L]
source(file.path(root, 'src/ldsc_gpca/r/gpsca_gwama_python_ldsc.r'))
fail <- function(expr, pattern) {
  error <- tryCatch({force(expr); NULL}, error = identity)
  stopifnot(inherits(error, 'error'), grepl(pattern, conditionMessage(error)))
}
traits <- c('C','A','B')
rows <- CJ(p1=traits,p2=traits)
rows[, `:=`(rg=ifelse(p1==p2,1,.4), se=.1, h2_obs=.2, h2_obs_se=.03,
             h2_int=1.05,h2_int_se=.01,gcov_int=ifelse(p1==p2,1.05,.1),gcov_int_se=.01)]
rows[, `:=`(z=rg/se,p=2*pnorm(-abs(rg/se)),rg_trait_wide=ifelse(p1==p2,1,.3),normalization_status='calculated')]
ldsc <- canonicalize_and_validate_ldsc(rows, traits)
select <- function(x=rows, method='trait_wide') select_correlation_matrix(ldsc,x,traits,method)
stopifnot(identical(select(method='pair'),ldsc$S_Stand), all(diag(select())==1),
          all(select()[upper.tri(select())]==.3), identical(rownames(select()),traits),
          identical(select(rows[.N:1]),select()), identical(select(rows[p1<=p2]),select()))
changed <- copy(rows); changed[, `:=`(se=.5,z=8,p=.8)]
stopifnot(identical(select(changed),select()),
          identical(compute_pc1(select(changed),traits)$loadings,compute_pc1(select(),traits)$loadings))
changed[1,rg_trait_wide:=NA_real_]
fail(select(changed),'unavailable')
stopifnot(identical(select(changed,'pair'),ldsc$S_Stand))
changed <- copy(rows); changed[1,normalization_status:='missing_reverse_pair']
fail(select(changed),'unavailable')
fail(select(rows[,!c('rg_trait_wide','normalization_status')]),'requires')
changed <- copy(rows); changed[p1=='A' & p2=='B',rg_trait_wide:=.5]
fail(select(changed),'Conflicting duplicate')
changed <- copy(rows); changed[p1==p2,rg_trait_wide:=.99]
fail(select(changed),'must equal 1')
changed <- copy(rows); changed[p1!=p2,rg_trait_wide:=1.1]
stopifnot(all(suppressWarnings(select(changed))[upper.tri(select())]==1.1))
fail(select_correlation_matrix(ldsc,changed,traits,'trait_wide',out_of_range_action='error'),'outside')
cat('Normalization selection, unavailable/duplicate/boundary values, order and default compatibility passed\n')
