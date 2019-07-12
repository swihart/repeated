# repeated R package
Bruce Swihart  
July 2019

## Re-Submission 1

  * https://github.com/swihart/repeated/issues/15
  * `repeated` was flagged for LTO mismatches using gcc9 and subsequently taken off CRAN.
  * LTOs are fixed with the following 3 edits:
  * 1.  In gausscop.f, line 157:  make 0 -> 0.0D0 for 6th argument in rs()
  * 2.  In chidden.f, line 436: comment out real precision rank and make rank integer in line above (435)
  * 3.  chidden.f:462:12:  `dqrcf(gmod,m,rank,qraux,work3,m,invec,info,1)` -> `dqrcf(gmod,m,rank,qraux,work3,m,invec,info)`
  
  In addition I've embellished the DESCRIPTION with information about models and citations to JK
  Lindsey's textbooks.

## Test environments
* local OS X install: R version 3.6.0 (2019-04-26)
* Ubuntu 14.04.5 LTS (on travis-ci): R version 3.6.0 (2017-01-27)
* win-builder: R Under development (unstable) (2019-07-05 r76784)

## R CMD check results
There were no ERRORs or WARNINGs or NOTEs.


## Downstream dependencies
There are currently no downstream dependencies for this package.

