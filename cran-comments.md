# repeated R package
Bruce Swihart  
January 2022


## Re-Submission 3

  * Corrected the LTO issues by correcting registration in repeated_init.c, as
  per BR email on 2022-01-31

## Re-Submission 2

  * https://github.com/swihart/repeated/issues/16
  * `call_R()` be gone!
  * rewrite the C-code of romberg integration in gnlmix to have entry .Call
  * `https:` where appropriate


## Test environments
* local OS X install: R version 4.1.2 (2021-11-01) 
    * Platform: x86_64-apple-darwin17.0 (64-bit)
    * Running under: macOS Mojave 10.14.6
* rhub: Apple Silicon (M1), macOS 11.6 Big Sur, R-release    
* rhub: Debian Linux, R-devel, clang, ISO-8859-15 locale
* rhub: Windows Server 2022, R-devel, 64 bit

## R CMD check results
There were no ERRORs or WARNINGs or NOTEs.


## Downstream dependencies
There are currently no downstream dependencies for this package.

