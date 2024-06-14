# repeated R package
Bruce Swihart  
Jun 2024

## Submission 1

  * fixed FORTRAN 2018 warnings in src/eigen.f as per Prof Ripley email
  * fixed FORTRAN portability: deprecated usage
  * removed logitord()
  
## Test environments
* local OS X install: R version 4.2.2 (2022-10-31)
    * Platform: x86_64-apple-darwin17.0 (64-bit)
    * Running under: macOS Big Sur 11.2.3
* rhub::check(platforms=c("gcc13","clang18"))
* 

## R CMD check results
There were no ERRORs or WARNINGs or NOTEs.


## Downstream dependencies
There are currently no downstream dependencies for this package.

# repeated R package
Bruce Swihart  
Mar 2022

## Submission

  * Corrected clang15 warnings and news.md formatting

## Test environments
* local OS X install: R version 4.1.2 (2021-11-01) 
    * Platform: x86_64-apple-darwin17.0 (64-bit)
    * Running under: macOS Mojave 10.14.6
* rhub: Debian Linux, R-devel, clang, ISO-8859-15 locale

## R CMD check results
There were no ERRORs or WARNINGs or NOTEs.


## Downstream dependencies
There are currently no downstream dependencies for this package.



# repeated R package
Bruce Swihart  
Mar 2022


## Submission

  * I messed up when I took over maintaining responsibilities -- the License is supposed to be GPL>=2 as Lindsey originally had.

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

