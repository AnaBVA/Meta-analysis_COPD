# https://seqqc.wordpress.com/2019/07/25/how-to-use-phyper-in-r/
# http://pedagogix-tagc.univ-mrs.fr/courses/ASG1/practicals/go_statistics_td/go_statistics_td_2015.html

## Compute the distribution of P-values (for all possible values of x).
x <- 15
m <- 1570
n <- 12929-1570
k <- 161
tail <- T # Under represented (T); Enrichment (F) 

p.value <- phyper(q=x-1, 
                   m=m, 
                   n=n, 
                   k=k, 
                   lower.tail= tail # Under represented (T); Enrichment (F) 
                   )

x.range <- phyper(1:(x*4), 
                  m=m, 
                  n=n, 
                  k=k, 
                  lower.tail= tail
     )

## Plot the P-value distribution on linear scales
plot(1:(x*4), 
     x.range, 
     type="l", 
     lwd=2, 
     col="violet", 
     main="Hypergeometric P-value for under represented num of genes", 
     xlab="Genes", 
     ylab="P-value = P[X > x]", 
     ylim=c(0, 1), 
     panel.first=grid()
     )

dens = dhyper(1:(x*4), m, n, k)

## Arrow indicating observed value
arrows(x, max(dens)*1.35, x, max(dens)*1.1, lwd=2, col='red', angle=30, length=0.1, code=2)
text(x, max(dens)*1.65, labels=paste("x=", x, "; p-val=", signif(digits=1, p.value), sep=""), col="red", font=2)

## We can plot the density below the P-value
lines(1:(x*4), dens, type="h", col="blue", lwd=2)
legend("topright", legend=c("P-value", "density"), lwd=2, col=c("violet", "blue"), bg="white", bty="o")




