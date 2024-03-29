library(rmutil)
library(repeated)

# Code to be run with
#   R -d "valgrind --tool=memcheck --leak-check=full --error-exitcode=1" --vanilla < tests/thisfile.R
# devtools::run_examples()

# Code to be run with
#   R -d "valgrind --tool=memcheck --leak-check=full --track-origins=yes" --vanilla < tests/thisfile.R

# dose <- c(9,12,4,9,11,10,2,11,12,9,9,9,4,9,11,9,14,7,9,8)
# y <- c(8.674419, 11.506066, 11.386742, 27.414532, 12.135699,  4.359469,
#        1.900681, 17.425948,  4.503345,  2.691792,  5.731100, 10.534971,
#       11.220260,  6.968932,  4.094357, 16.393806, 14.656584,  8.786133,
#       20.972267, 17.178012)
# resp <- restovec(matrix(y, nrow=4, byrow=TRUE), name="y")
# reps <- rmna(resp, tvcov=tvctomat(matrix(dose, nrow=4, byrow=TRUE), name="dose"))

# same linear normal model with random normal intercept fitted four ways
# compare with growth::elliptic(reps, model=~dose, preg=c(0,0.6), pre=4)
# glmm(y~dose, nest=individuals, data=reps)
# gnlmm(reps, mu=~dose, pmu=c(8.7,0.25), psh=3.5, psd=3)
#gnlmix(reps, mu=~a+b*dose+rand, random="rand", pmu=c(8.7,0.25),
#            pshape=3.44, pmix=2.3)


library(repeated)

dose <- c(9,12,4,9,11,10,2,11,12,9,9,9,4,9,11,9,14,7,9,8)
y <- c(8.674419, 11.506066, 11.386742, 27.414532, 12.135699,  4.359469,
       1.900681, 17.425948,  4.503345,  2.691792,  5.731100, 10.534971,
       11.220260,  6.968932,  4.094357, 16.393806, 14.656584,  8.786133,
       20.972267, 17.178012)
id <- rep(1:4, each=5)


beg.ex1 <- Sys.time()
ex1 <- gnlmix(y, mu=~a+b*dose+rand, random="rand", nest=id, pmu=c(a=8.7,b=0.25),
              pshape=3.44, pmix=2.3)
end.ex1 <- Sys.time()
time.ex1 <- end.ex1 - beg.ex1

ex1$coef

ex1$maxlike

time.ex1
