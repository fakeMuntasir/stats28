library(tidyverse)

#' Test of hypothesis which determines the p-value,degrees of freedom,confidence interval and test
#' statistic given the arguments
#'
#' @param approach argument to determine whether the test is parametric or non-parametric
#' @param independent For matched pair data set, this parameter is false, and otherwise true.
#' @param x A random variable from a normal distribution
#' @param y Another random variable from a distribution
#' @param mu_0 mean of a specific distribution
#' @param sigma_x population standard deviation of given data x
#' @param sigma_y population standard deviation of given data y
#' @param sigma_d population standard deviation of paired differences
#' @param n1  number of rows matrix
#' @param n2 number of columns matrix
#' @param B number of replication
#' @param alpha level of significance
#' @param H1 alternate hypothesis
#'
#' @return p-value,degrees-of-freedom,confidence-interval,test-statistic
#'
#' @examples test_of_hypothesis(approach = "parametric", independent = FALSE, x = x,
#' y = y, mu_0 = mu_0, alpha = alpha, H1 = H1 )



test_of_hypothesis<-function(approach,independent=NULL,x,y=NULL,mu_0=NULL,sigma_x=NULL,
                             sigma_y=NULL,sigma_d=NULL,n1=NULL,n2=NULL,B=NULL,alpha,H1){

  if(approach == "parametric"){



    if( is.vector(x) & !is.character(x) & !is.null(mu_0)){

      mu=mean(x)
      n=length(x)


      if(!is.null(sigma_x) & n>30){

        z_cal=round(((mu-mu_0)*sqrt(n))/sigma_x,3)
        z_tab<-qnorm(alpha,0,1,F)
        z_tab2<-qnorm(alpha/2,0,1,F)


        b<-function(H1){
          if(H1=="mu>mu_0"){
            p=round(pnorm(z_cal,lower.tail = F),3)
            lo=NULL
            up=round(mu+((z_tab*sigma)/sqrt(n)),3)

          }
          else if(H1=="mu<mu_0"){
            p=round(pnorm(z_cal,lower.tail = T),3)
            lo=round(mu-((z_tab*sigma)/sqrt(n)),3)
            up=NULL


          }
          else if(H1== "mu!=mu_0"){
            p=round(2*(1-pnorm(abs(z_cal))),3)
            lo=round(mu-((z_tab2*sigma)/sqrt(n)),3)
            up=round(mu+((z_tab2*sigma)/sqrt(n)),3)


          }
          H0<-"mu=mu_0"


          return(tibble(p_value=p,Test_statistics=z_cal,lower_bound=lo,upper_bound=up ,H1=H1))
        }
      }
      else if(n<30){
        s=sqrt((sum((x-mu)^2))/(n-1))

        t_cal=round(((mu-mu_0)*sqrt(n))/s,3)
        t_tab<-qt((alpha),(n-1),lower.tail = F)
        t_tab_2<-qt((alpha/2),(n-1),lower.tail = F)

        b<-function(H1){
          if(H1=="mu>mu_0"){
            p=round(pt(t_cal,(n-1),lower.tail = F),3)
            lo=NULL
            up=round(mu+((t_tab*s)/sqrt(n)),3)


          }
          else if(H1=="mu<mu_0"){
            p=round(pt(t_cal,(n-1),lower.tail = T),3)
            lo=round(mu-((t_tab*s)/sqrt(n)),3)
            up=NULL


          }
          else if(H1== "mu!=mu_0"){
            p=round(2*pt(-abs(t_cal),(n-1)),3)
            lo=round(mu-((t_tab_2*s)/sqrt(n)),3)
            up=round(mu+((t_tab_2*s)/sqrt(n)),3)


          }
          H0<-"mu=mu_0"
          df<-n-1
          return(tibble(p_value=p,Test_statistics=t_cal,degrees_of_freedom=df,lower_bound=lo,upper_bound=up,H1=H1))
        }
      }
      return(b(H1))

    }


    if(!is.null(independent)  &  !is.null(y)){

      n=length(x)
      m=length(y)


      if(independent == "TRUE" ){

        if( !is.null(sigma_x) & !is.null(sigma_y) & n>30 & m>30){

          xbar <- mean(x)
          ybar <- mean(y)

          co<-sqrt(((sigma_x)^2/n)+((sigma_y)^2/m))
          z_cal<-round((xbar-ybar)/co,3)

          z_tab<-qnorm(alpha,0,1,F)
          z_tab2<-qnorm(alpha/2,0,1,F)


          b<-function(H1){
            if(H1=="mu>mu_0"){
              p=round(pnorm(z_cal,lower.tail = F),3)
              lo=NULL
              up=round((xbar-ybar)+(z_tab*co),3)

            }
            else if(H1=="mu<mu_0"){
              p=round(pnorm(z_cal,lower.tail = T),3)
              lo=round((xbar-ybar)-(z_tab*co),3)
              up=NULL


            }
            else if(H1== "mu!=mu_0"){
              p=round(2*(1-pnorm(abs(z_cal))),3)
              lo=round((xbar-ybar)-(z_tab2*co),3)
              up=round((xbar-ybar)+(z_tab2*co),3)


            }
            H0<-"mu=mu_0"


            return(tibble(p_value=p,Test_statistics=z_cal,lower_bound=lo,upper_bound=up,H1=H1))
          }


        }
        if( is.null(sigma_x) & is.null(sigma_y) & n<30 & m<30){

          xbar <- mean(x)
          ybar <- mean(y)

          sx <- sd(x)
          sy <- sd(y)
          sp <- sqrt(((n-1)*sx^2+(m-1)*sy^2)/(n+m-2))
          t_cal <- round((xbar - ybar) / sp*(sqrt(1/n+1/m)),3)
          t_tab <-qt(alpha, (n+m-2),lower.tail = F)
          t_tab2 <- qt((alpha/2), (n+m-2),lower.tail = F)


          b<-function(H1){
            if(H1=="mu>mu_0"){
              p=round(pt(t_cal,(n+m-2),lower.tail = F),3)
              lo=NULL
              up=round((xbar-ybar)+(t_tab*sp*(sqrt(1/n+1/m))),3)


            }
            else if(H1=="mu<mu_0"){
              p=round(pt(t_cal,(n+m-2),lower.tail = T),3)
              lo=round((xbar-ybar)-(t_tab*sp*(sqrt(1/n+1/m))),3)
              up=NULL


            }
            else if(H1== "mu!=mu_0"){
              p=round(2*pt(-abs(t_cal),(n+m-2)),3)
              lo=round((xbar-ybar)-(t_tab2*sp*(sqrt(1/n+1/m))),3)
              up=round((xbar-ybar)+(t_tab2*sp*(sqrt(1/n+1/m))),3)


            }
            H0<-"mu=mu_0"
            df<-n+m-2
            return(tibble(p_value=p,Test_statistics=t_cal,degrees_of_freedom=df,lower_bound=lo,upper_bound=up,H1=H1))
          }
        }

      }

      else if(independent == "FALSE")
      {

        if( !is.null(sigma_d) & n>30 & m>30){

          differences <- (y -x)
          r <- length(differences)

          mean_diff <- mean(differences)

          co<-sqrt(((sigma_d)^2)/r)
          z_cal<-round(mean_diff/co,3)

          z_tab<-qnorm(alpha,0,1,F)
          z_tab2<-qnorm((alpha/2),0,1,F)


          b<-function(H1){
            if(H1=="mu>mu_0"){
              p=round(pnorm(z_cal,lower.tail = F),3)
              lo=NULL
              up=round(mean_diff+((z_tab*sigma_d)/sqrt(r)),3)

            }
            else if(H1=="mu<mu_0"){
              p=round(pnorm(z_cal,lower.tail = T),3)
              lo=round(mean_diff-((z_tab*sigma_d)/sqrt(r)),3)
              up=NULL


            }
            else if(H1== "mu!=mu_0"){
              p=round(2*(1-pnorm(abs(z_cal))),3)
              lo=round(mean_diff-((z_tab2*sigma_d)/sqrt(r)),3)
              up=round(mean_diff+((z_tab*sigma_d)/sqrt(r)),3)


            }
            H0<-"mu=mu_0"


            return(tibble(p_value=p,Test_statistics=z_cal,lower_bound=lo,upper_bound=up ,H1=H1))
          }


        }

        if( is.null(sigma_d) & n<30 & m<30){


          differences <- (y -x)
          r <- length(differences)

          mean_diff <- mean(differences)
          sd_diff <- sqrt(sum((differences - mean_diff)^2) / (r - 1))

          t_cal <- round(mean_diff / (sd_diff / sqrt(r)),3)

          df <- r- 1

          t_tab <-qt(alpha, df,lower.tail = F)
          t_tab2 <- qt((alpha/2),df,lower.tail = F)


          b<-function(H1){
            if(H1=="mu>mu_0"){
              p=round(pt(t_cal,df,lower.tail = F),3)
              lo=NULL
              up=round(mean_diff+((t_tab*sd_diff)/sqrt(r)),3)


            }
            else if(H1=="mu<mu_0"){
              p=round(pt(t_cal,df,lower.tail = T),3)
              lo=round(mean_diff-((t_tab*sd_diff)/sqrt(r)),3)
              up=NULL


            }
            else if(H1== "mu!=mu_0"){
              p=round(2*pt(-abs(t_cal),df),3)
              lo=round(mean_diff-((t_tab2*sd_diff)/sqrt(r)),3)
              up=round(mean_diff+((t_tab2*sd_diff)/sqrt(r)),3)


            }
            H0<-"mu=mu_0"

            return(tibble(p_value=p,Test_statistics=t_cal,degrees_of_freedom=df,lower_bound=lo,upper_bound=up,H1=H1))
          }

        }

      }

      return(b(H1))

    }


    if(is.character(x) & !is.null(y) & !is.null(n1) & !is.null(n2)){



      datt <- data.frame(x,y )

      N <- n1*n2

      constant <- sum(y)^2 / N

      sst <- sum(y^2) - constant

      yi. <- sapply(split(y,x), sum )
      sstreat <- (sum(yi.^2)/n2) - constant

      sse <- sst - sstreat

      mstreat <- sstreat / (n1-1)
      mse <- sse / (N-n1)

      f_cal <- round(mstreat / mse,3)

      f_tab <- qf((alpha), n1-1, N-n1, lower.tail = F)



      H0 <- "Treatment means are equal"

      pvalue <-round(pf(f_cal,n1-1, N-n1,lower.tail = F),3)
      t_tab<-qt((alpha/2),(N-n1),lower.tail = F)
      t_tab
      val<-t_tab*(sqrt(mse/n2))
      val
      a<-datt %>%
        group_by(x) %>%
        reframe(CI_lower=(mean(y)-val),CI_upper=(mean(y)+val))
      df1<-n1-1
      df2<-N-n1



      return(data.frame(p_value=pvalue,test_statistic =  f_cal,df1=df1,df2=df2,a,H1=H1))
    }


    if(is.matrix(x) & !is.null(n1)  & !is.null(n2)){



      df=(n1-1)*(n2-1)



      a<-NULL
      b<-NULL
      for(i in 1:n1){
        a[i]<-sum(x[i,])
      }

      for(j in 1:n2){
        b[j]<-sum(x[,j])
      }


      eij <- matrix(0, nrow = n1, ncol = n2)
      for(i in 1:n1){
        for(j in 1:n2){
          eij[i,j]<-a[i]*b[j]/sum(x)
        }
      }

      q<-0

      for(i in 1:n1){
        for(j in 1:n2){
          q=q+sum((x[i,j]-eij[i,j])^2/eij[i,j])
        }
      }
      r<-round(q,3)
      pvalue<-round(pchisq(q,df,lower.tail = F),3)

      H0<-"No association"


      return(data.frame(p_value=pvalue,Test_statistic=r,degrees_of_freedom=df,H1=H1))

    }
  }



  if(approach == "non_parametric"){

    if( !is.null(y) & !is.null(B)){

      boot_stat_vec <- c()

      for(i in 1:B){

        dat <- data.frame(y, x)

        index <- sample(1:length(x), size = length(x), replace = T)

        dat_boot <- dat[index,]

        boot_stat_vec[i] <-  dat_boot %>%
          group_by(x) %>%
          summarise(mean = mean(y)) %>%
          summarise(test_stat = mean[x == 1] - mean[x == 0]) %>%
          pull()
      }
      test_stat_main <- dat %>%
        group_by(x) %>%
        summarise(mean = mean(y)) %>%
        summarise(test_stat = mean[x == 1] - mean[x == 0]) %>%
        pull()

      p_val = (sum(boot_stat_vec >= test_stat_main))/B

      H0<-"There is no treatment effect"

      return(data.frame(p_val=p_val,test_statistic=test_stat_main,H1=H1))

    }
  }

}

