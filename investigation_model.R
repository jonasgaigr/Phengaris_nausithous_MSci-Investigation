library(tidyverse)
library(sf)
library(sfnetworks)
library(units)
library(tidygraph)
library(leaflet)
library(sp)
library(rjags)
library(runjags)

# LOAD DATA
# Dataset of patches
patch_sjtsk <- sf::st_read("https://github.com/jonasgaigr/Phengaris_nausithous_MSci-Investigation/blob/main/patches_sjtsk_finale.gpkg?raw=true")
# Study area polygon
krivoklatsko_sjtsk <- st_read("https://github.com/jonasgaigr/Phengaris_nausithous_MSci-Investigation/blob/main/krivoklatsko_sjtsk.gpkg?raw=true")

# Define priors ----
tmin <- 25
tmax <- tmin+25
a0 <- -0.5 # Baseline persistence
a1 <- 1 # Effect of host plant on persistance
a2 <- 0.05 # Effect of patch size on persistance
b0 <- 1 # Baseline colonisation
b1 <- -0.005 # Effect of distance on colonisation

# Metapopulation model ----
no <- nrow(patch_sjtsk)
s <- matrix(0, no, tmax)
s[,1] <- as.numeric(patch_sjtsk$INIT_OCC)
A <- patch_sjtsk %>% dplyr::pull(SHAPEAREA)/10000
P <- patch_sjtsk %>% dplyr::pull(PLANT)

for(t in 1:(tmax-1)) {
  for(i in 1:no) {
    if(s[i,t] == 1) {
      # Persistance
      l <- a0+a1*P[i]+a2*A[i]
      q <- exp(l)/(1+exp(l))
      s[i,t+1] <- rbinom(1,1,q)
    } else {
      # Colonisation
      pj <- rep(0,no)
      for(j in 1:no) {
        l <- b0+b1*dist[j,i]
        pj[j] <- 1-exp(l)/(1+exp(l))*s[j,t]
      }
      p <- 1-prod(pj)*0.999999
      s[i, t+1] <- rbinom(1,1,p)
    }
  }
}

Oc <- s[,tmax]
si <- s*NA
si[,tmax] <- as.numeric(patch_sjtsk$INIT_OCC)

# JAGS Code ----
butterflies <- "model

{ 

for(t in 1:(tmax-1)) {
  for(i in 1:no) {
  # Persistence
    logit(q[i,t]) <- a0+a1*P[i]+a2*A[i]
    # Colonisation
    for(j in 1:no) {
      l[t,i,j] <- b0-b1*dist[j,i]
      pj[t,i,j] <- 1-(exp(l[t,i,j])/(1+exp(l[t,i,j]))*si[j,t]) # Failure to colonise
    }
  p[i,t] <- 1-prod(pj[t,i,1:no])*0.999999
  r[i,t] <-ifelse(si[i,t] == 1, q[i,t], p[i,t])
  si[i,t+1]~dbin(r[i,t],1)
  }
}

# Occupancy data
for(i in 1:no) {
    poc[i]<-sum(si[i,tmin:tmax])/(tmax-tmin+1)+0.00001
    Oc[i]~dbin(poc[i],1)
}
  #data# tmin,tmax,no,P,A,dist,Oc
  #monitor# a0,a1,a2,b0,b1
  
  #### PRIORS ####
  a0~dnorm(0,1)
  a1~dnorm(0,1)
  a2~dnorm(0,1)
  
  b0~dnorm(0,1)
  b1~dgamma(0.005,0.005)
  
  for(i in 1:no) {si[i,1]~dbin(0.5,1)}
}

"

# Explore the change of precision with increasing sample size (plot median, L95, u95) ----
patch_6 <- patch_sjtsk %>%
  st_filter(., st_buffer(st_centroid(krivoklatsko_sjtsk), dist = 3300))
patch_12 <- patch_sjtsk %>%
  st_filter(., st_buffer(st_centroid(krivoklatsko_sjtsk), dist = 4200))
patch_24 <- patch_sjtsk %>%
  st_filter(., st_buffer(st_centroid(krivoklatsko_sjtsk), dist = 5400))
patch_31 <- patch_sjtsk %>%
  st_filter(., st_buffer(st_centroid(krivoklatsko_sjtsk), dist = 6000))
patch_40 <- patch_sjtsk %>%
  st_filter(., st_buffer(st_centroid(krivoklatsko_sjtsk), dist = 7000))
patch_50 <- patch_sjtsk %>%
  st_filter(., st_buffer(st_centroid(krivoklatsko_sjtsk), dist = 9000))
patch_65 <- patch_sjtsk %>%
  st_filter(., st_buffer(st_centroid(krivoklatsko_sjtsk), dist = 11000))
patch_82 <- patch_sjtsk %>%
  st_filter(., st_buffer(st_centroid(krivoklatsko_sjtsk), dist = 13500))


# 6 patches ----
dist <- sf::st_distance(patch_6) %>%
  units::drop_units() 

no <- nrow(patch_6)
s <- matrix(0, no, tmax)
s[,1] <- as.numeric(patch_6$INIT_OCC)
s[,1]
s[,1] <- rbinom(no, 1, 0.5)
A <- patch_6 %>% dplyr::pull(SHAPEAREA)/10000
P <- patch_6 %>% dplyr::pull(PLANT)

# Main simulation
for(t in 1:(tmax-1)) {
  for(i in 1:no) {
    if(s[i,t] == 1) {
      # Persistance
      l <- a0+a1*P[i]+a2*A[i]
      q <- exp(l)/(1+exp(l))
      s[i,t+1] <- rbinom(1,1,q)
    } else {
      # Colonisation
      pj <- rep(0,no)
      for(j in 1:no) {
        l <- b0+b1*dist[j,i]
        pj[j] <- 1-exp(l)/(1+exp(l))*s[j,t]
      }
      p <- 1-prod(pj)*0.999999
      s[i, t+1] <- rbinom(1,1,p)
    }
  }
}

Oc <- s[,tmax]
si <- s*NA
si[,tmax] <- s[,tmax]
si[,tmax] <- as.numeric(patch_24$INIT_OCC)

results_6 <- run.jags(butterflies, burnin = 500,
                    sample = 1000,
                    adapt = 500,
                    n.chains = 2)
sumr_6 <- summary(results_6)
sumr_6
plot(results_6, c("trace","hist"))

a0 <- sumr_6[1,1]
a1 <- sumr_6[2,1]
a2 <- sumr_6[3,1]
b0 <- sumr_6[4,1]
b1 <- sumr_6[5,1]


# 12 patches ----
dist <- sf::st_distance(patch_12) %>%
  units::drop_units() 
dist <- dist+0.0000001
no <- nrow(patch_12)
s <- matrix(0, no, tmax)
s[,1] <- as.numeric(patch_12$INIT_OCC)
#s[,1] <- rbinom(no, 1, 0.5)
A <- patch_12 %>% dplyr::pull(SHAPEAREA)/10000
P <- patch_12 %>% dplyr::pull(PLANT)
s[,1]
#P <- rbinom(no, 1, 0.5)

# Main simulation
for(t in 1:(tmax-1)) {
  for(i in 1:no) {
    if(s[i,t] == 1) {
      # Persistance
      l <- a0+a1*P[i]+a2*A[i]
      q <- exp(l)/(1+exp(l))
      s[i,t+1] <- rbinom(1,1,q)
    } else {
      # Colonisation
      pj <- rep(0,no)
      for(j in 1:no) {
        l <- b0+b1*dist[j,i]
        pj[j] <- 1-exp(l)/(1+exp(l))*s[j,t]
      }
      p <- 1-prod(pj)*0.999999
      s[i, t+1] <- rbinom(1,1,p)
    }
  }
}

a0 <- sumr_12[1,1]
a1 <- sumr_12[2,1]
a2 <- sumr_12[3,1]
b0 <- sumr_12[4,1]
b1 <- sumr_12[5,1]

Oc <- s[,tmax]
si <- s*NA
si[,tmax] <- s[,tmax]
si[,tmax] <- as.numeric(patch_12$INIT_OCC)
si[,tmax]
s

results_12 <- run.jags(butterflies, burnin = 500,
                       sample = 1000,
                       adapt = 500,
                       n.chains = 2)
sumr_12 <- summary(results_12)
sumr_12
plot(results_12, c("trace","hist"))
drop.k(results_12, drop = "A[1:5]", k=1)


# 24 patches ----
dist <- sf::st_distance(patch_24) %>%
  units::drop_units() 
dist <- dist
no <- nrow(patch_24)
s <- matrix(0, no, tmax)
s[,1] <- as.numeric(patch_24$INIT_OCC)
s[,1]
#s[,1] <- rbinom(no, 1, 0.5)
A <- patch_24 %>% dplyr::pull(SHAPEAREA)/10000
P <- patch_24 %>% dplyr::pull(PLANT)
P

# Main simulation
for(t in 1:(tmax-1)) {
  for(i in 1:no) {
    if(s[i,t] == 1) {
      # Persistance
      l <- a0+a1*P[i]+a2*A[i]
      q <- exp(l)/(1+exp(l))
      s[i,t+1] <- rbinom(1,1,q)
    } else {
      # Colonisation
      pj <- rep(0, no)
      for(j in 1:no) {
        l <- b0 + b1*dist[j, i]
        pj[j] <- 1-exp(l)/(1+exp(l))*s[j,t]
      }
      p <- 1 - prod(pj)*0.999999
      s[i, t+1] <- rbinom(1,1,p)
    }
  }
}

Oc <- s[,tmax]
si <- s*NA
si[,tmax] <- s[,tmax]
si[,1] <- as.numeric(patch_24$INIT_OCC)
#si[,1] <- s[,tmax]
si
s

results_24 <- run.jags(butterflies, burnin = 500,
                       sample = 1000,
                       adapt = 500,
                       n.chains = 2)
sumr_24 <- summary(results_24)
sumr_24
plot(results_24, c("trace","hist"))

# 31 patches ----
dist <- sf::st_distance(patch_31) %>%
  units::drop_units() 
no <- nrow(patch_31)
s <- matrix(0, no, tmax)
s[,1] <- as.numeric(patch_31$INIT_OCC)
s[,1]
#s[,1] <- rbinom(no, 1, 0.5)
A <- patch_31 %>% dplyr::pull(SHAPEAREA)/10000
P <- patch_31 %>% dplyr::pull(PLANT)

# Main simulation
for(t in 1:(tmax-1)) {
  for(i in 1:no) {
    if(s[i,t] == 1) {
      # Persistance
      l <- a0+a1*P[i]+a2*A[i]
      q <- exp(l)/(1+exp(l))
      s[i,t+1] <- rbinom(1,1,q)
    } else {
      # Colonisation
      pj <- rep(0,no)
      for(j in 1:no) {
        l <- b0+b1*dist[j,i]
        pj[j] <- 1-exp(l)/(1+exp(l))*s[j,t]
      }
      p <- 1-prod(pj)*0.999999
      s[i,t+1] <- rbinom(1,1,p)
    }
  }
}

rowSums(s)/tmax
Oc <- s[,tmax]
si <- s*NA
#si[,tmax] <- s[,tmax]
si[,tmax] <- as.numeric(patch_31$INIT_OCC)
s

results_31 <- run.jags(butterflies, burnin = 500,
                       sample = 1000,
                       adapt = 500,
                       n.chains = 2)
sumr_31 <- summary(results_31)
sumr_31
plot(results_31, c("trace","hist"))


# 40 patches ----
dist <- sf::st_distance(patch_40) %>%
  units::drop_units() 
dist
no <- nrow(patch_40)
s <- matrix(0, no, tmax)
s[,1] <- as.numeric(patch_40$INIT_OCC)
#s[,1] <- rbinom(no, 1, 0.5)
A <- patch_40 %>% dplyr::pull(SHAPEAREA)/10000
P <- patch_40 %>% dplyr::pull(PLANT)
#P <- rbinom(no, 1, 0.5)

# Main simulation
for(t in 1:(tmax-1)) {
  for(i in 1:no) {
    if(s[i,t] == 1) {
      # Persistance
      l <- a0+a1*P[i]+a2*A[i]
      q <- exp(l)/(1+exp(l))
      s[i,t+1] <- rbinom(1, 1, q)
    } else {
      # Colonisation
      pj <- rep(0, no)
      for(j in 1:no) {
        l <- b0+b1*dist[j, i]
        pj[j] <- 1-exp(l)/(1+exp(l))*s[j,t]
      }
      p <- 1 - prod(pj)*0.999999
      s[i, t+1] <- rbinom(1, 1, p)
    }
  }
}

rowSums(s)/tmax
Oc <- s[,tmax]
si <- s*NA
#si[,tmax] <- s[,tmax]
si[,tmax] <- as.numeric(patch_40$INIT_OCC)
s
si

results_40 <- run.jags(butterflies, burnin = 500,
                       sample = 1000,
                       adapt = 500,
                       n.chains = 2)
sumr_40 <- summary(results_40)
sumr_40
plot(results_40, c("trace","hist"))



# 50 patches ----
dist <- sf::st_distance(patch_50) %>%
  units::drop_units() 
no <- nrow(patch_50)
s <- matrix(0, no, tmax)
s[,1] <- as.numeric(patch_50$INIT_OCC)
#s[,1] <- rbinom(no, 1, 0.5)
A <- patch_50 %>% dplyr::pull(SHAPEAREA)/10000
P <- patch_50 %>% dplyr::pull(PLANT)
A %>% hist(., breaks = 100)
# Main simulation
for(t in 1:(tmax-1)) {
  for(i in 1:no) {
    if(s[i,t] == 1) {
      # Persistance
      l <- a0+a1*P[i]+a2*A[i]
      q <- exp(l)/(1+exp(l))
      s[i,t+1] <- rbinom(1,1,q)
    } else {
      # Colonisation
      pj <- rep(0,no)
      for(j in 1:no) {
        l <- b0+b1*dist[j,i]
        pj[j] <- 1-exp(l)/(1+exp(l))*s[j,t]
      }
      p <- 1 - prod(pj)*0.999999
      s[i,t+1] <- rbinom(1,1,p)
    }
  }
}

Oc <- s[,tmax]
si <- s*NA
#si[,tmax] <- s[,tmax]
si[,tmax] <- as.numeric(patch_50$INIT_OCC)
si
s

results_50 <- run.jags(butterflies, burnin = 500,
                       sample = 1000,
                       adapt = 500,
                       n.chains = 2)
sumr_50 <- summary(results_50)
sumr_50
plot(results_50, c("trace","hist"))

# 65 patches ----
dist <- sf::st_distance(patch_65) %>%
  units::drop_units() 
no <- nrow(patch_65)
s <- matrix(0, no, tmax)
s[,1] <- as.numeric(patch_65$INIT_OCC)
s[,1] <- rbinom(no, 1, 0.5)
A <- patch_65 %>% dplyr::pull(SHAPEAREA)/10000
P <- patch_65 %>% dplyr::pull(PLANT)

# Main simulation
for(t in 1:(tmax-1)) {
  for(i in 1:no) {
    if(s[i,t] == 1) {
      # Persistance
      l <- a0 + a1*P[i] + a2*A[i]
      q <- exp(l)/(1+exp(l))
      s[i,t+1] <- rbinom(1, 1, q)
    } else {
      # Colonisation
      pj <- rep(0, no)
      for(j in 1:no) {
        l <- b0 + b1*dist[j, i]
        pj[j] <- 1-exp(l)/(1+exp(l))*s[j,t]
      }
      p <- 1 - prod(pj)*0.999999
      s[i, t+1] <- rbinom(1, 1, p)
    }
  }
}

Oc <- s[,tmax]
si <- s*NA
si[,tmax] <- s[,tmax]

results_65 <- run.jags(butterflies, burnin = 500,
                       sample = 1000,
                       adapt = 500,
                       n.chains = 2)
sumr_65 <- summary(results_65)


# 82 patches ----
dist <- sf::st_distance(patch_82) %>%
  units::drop_units() 
patch_82
no <- nrow(patch_82)
s <- matrix(0, no, tmax)
as.numeric(patch_82$INIT_OCC)
s[,1] <- as.numeric(patch_40$INIT_OCC)
#s[,1] <- rbinom(no, 1, 0.5)
A <- patch_82 %>% dplyr::pull(SHAPEAREA)/10000
P <- patch_82 %>% dplyr::pull(PLANT)

# Main simulation
for(t in 1:(tmax-1)) {
  for(i in 1:no) {
    if(s[i,t] == 1) {
      # Persistance
      l <- a0 + a1*P[i] + a2*A[i]
      q <- exp(l)/(1+exp(l))
      s[i,t+1] <- rbinom(1, 1, q)
    } else {
      # Colonisation
      pj <- rep(0, no)
      for(j in 1:no) {
        l <- b0 + b1*dist[j, i]
        pj[j] <- 1-exp(l)/(1+exp(l))*s[j,t]
      }
      p <- 1 - prod(pj)*0.999999
      s[i, t+1] <- rbinom(1, 1, p)
    }
  }
}

Oc <- s[,tmax]
si <- s*NA
si[,tmax] <- s[,tmax]

results_82 <- run.jags(butterflies, burnin = 500,
                       sample = 1000,
                       adapt = 500,
                       n.chains = 2)
sumr_82 <- summary(results_82)

# ALL patches ----
patch_sjtsk
dist <- sf::st_distance(patch_sjtsk) %>%
  units::drop_units() 
dist <- dist
dist
no <- nrow(patch_sjtsk)
s <- matrix(0, no, tmax)
s[,1] <- as.numeric(patch_sjtsk$INIT_OCC)
A <- patch_sjtsk %>% dplyr::pull(SHAPEAREA)/10000
P <- patch_sjtsk %>% dplyr::pull(PLANT)

# Main simulation
for(t in 1:(tmax-1)) {
  for(i in 1:no) {
    if(s[i,t] == 1) {
      # Persistance
      l <- a0+a1*P[i]+a2*A[i]
      q <- exp(l)/(1+exp(l))
      s[i,t+1] <- rbinom(1,1,q)
    } else {
      # Colonisation
      pj <- rep(0,no)
      for(j in 1:no) {
        l <- b0+b1*dist[j,i]
        pj[j] <- 1-exp(l)/(1+exp(l))*s[j,t]
      }
      p <- 1-prod(pj)*0.999999
      s[i, t+1] <- rbinom(1,1,p)
    }
  }
}

rowSums(s)/tmax
Oc <- s[,tmax]
si <- s*NA
#si[,tmax] <- s[,tmax]
si[,tmax] <- as.numeric(patch_sjtsk$INIT_OCC)
si[,tmax]
s

results <- run.jags(butterflies, burnin = 500,
                    sample = 1000,
                    adapt = 500,
                    n.chains = 2)
sumr <- summary(results)
plot(results, c("trace","hist"))


# Plot precision change ----

sumr

sumr_overall <- bind_rows(sumr_6 %>%
                            as.data.frame() %>%
                            mutate(patches = 6,
                                   parameter = c("a0", "a1", "a2", "b0", "b1")),
                          sumr_12 %>%
                            as.data.frame() %>%
                            mutate(patches = 12,
                                   parameter = c("a0", "a1", "a2", "b0", "b1")),
                          sumr_24 %>%
                            as.data.frame() %>%
                            mutate(patches = 24,
                                   parameter = c("a0", "a1", "a2", "b0", "b1")),
                          sumr_31 %>%
                            as.data.frame() %>%
                            mutate(patches = 31,
                                   parameter = c("a0", "a1", "a2", "b0", "b1")),
                          sumr_40 %>%
                            as.data.frame() %>%
                            mutate(patches = 40,
                                   parameter = c("a0", "a1", "a2", "b0", "b1")),
                          sumr_50 %>%
                            as.data.frame() %>%
                            mutate(patches = 50,
                                   parameter = c("a0", "a1", "a2", "b0", "b1")),
                          sumr_65 %>%
                            as.data.frame() %>%
                            mutate(patches = 65,
                                   parameter = c("a0", "a1", "a2", "b0", "b1")),
                          sumr_82 %>%
                            as.data.frame() %>%
                            mutate(patches = 82,
                                   parameter = c("a0", "a1", "a2", "b0", "b1")),
                          sumr %>%
                            as.data.frame() %>%
                            mutate(patches = 90,
                                   parameter = c("a0", "a1", "a2", "b0", "b1")
                                   ))

sumr_overall

ggplot(data = sumr_overall, aes(x = patches, 
                                y = Median)) +
  geom_point() +
  geom_errorbar(aes(ymin = Lower95, ymax = Upper95)) +
  facet_wrap(~parameter) +
  labs(xlab = "number pf patcehs",
       caption = "a0 = -0.2; a1 = 0.5; a2 = 1; b0 = 0; b1 = -0.001") +
  theme_bw()


plot_a0 <- ggplot(data = filter(sumr_overall, parameter == "a0"), aes(x = patches, 
                                y = Median)) +
  geom_point() +
  geom_errorbar(aes(ymin = Lower95, ymax = Upper95)) +
  geom_hline((aes(yintercept = a0))) +
  xlab("number pf patches") +
  ggtitle("a0") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) 
plot_a0
plot_a1 <- ggplot(data = filter(sumr_overall, parameter == "a1"), aes(x = patches, 
                                                                   y = Median)) +
  geom_point() +
  geom_errorbar(aes(ymin = Lower95, ymax = Upper95)) +
  geom_hline((aes(yintercept = a1))) +
  xlab("number pf patches") +
  ggtitle("a1") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
plot_a1
plot_a2 <- ggplot(data = filter(sumr_overall, parameter == "a2"), aes(x = patches, 
                                                                   y = Median)) +
  geom_point() +
  geom_errorbar(aes(ymin = Lower95, ymax = Upper95)) +
  geom_hline((aes(yintercept = a2))) +
  xlab("number pf patches") +
  ggtitle("a2") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
plot_a2
plot_b0 <- ggplot(data = filter(sumr_overall, parameter == "b0"), aes(x = patches, 
                                                                   y = Median)) +
  geom_point() +
  geom_errorbar(aes(ymin = Lower95, ymax = Upper95)) +
  geom_hline((aes(yintercept = b0))) +
  xlab("number pf patches") +
  ggtitle("b0") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
plot_b0
plot_b1 <- ggplot(data = filter(sumr_overall, parameter == "b1"), aes(x = patches, 
                                                                   y = Median)) +
  geom_point() +
  geom_errorbar(aes(ymin = Lower95, ymax = Upper95)) +
  geom_hline((aes(yintercept = b1))) +
  scale_y_continuous(breaks = c(3.5,2.5), limits = c(0,3)) +
  xlab("number pf patches") +
  ggtitle("b1") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
plot_b1

#Plot occuapancy pattern
patch_sjtsk <- patch_sjtsk %>%
  mutate(colour = case_when(INIT_OCC == 1 ~ "blue",
                            INIT_OCC == 0 ~ "red"))

bbox_new <- st_bbox(krivoklatsko_sjtsk) # current bounding box
xrange <- bbox_new$xmax - bbox_new$xmin # range of x values
yrange <- bbox_new$ymax - bbox_new$ymin # range of y values
bbox_new[1] <- bbox_new[1] - (0.025 * xrange) # xmin - left
bbox_new[3] <- bbox_new[3] + (0.025 * xrange) # xmax - right
bbox_new[2] <- bbox_new[2] - (0.025 * yrange) # ymin - bottom
bbox_new[4] <- bbox_new[4] + (0.025 * yrange) # ymax - top

bbox_new <- bbox_new %>%  # take the bounding box ...
  st_as_sfc()
study_map <- ggplot2::ggplot() +
  ggplot2::geom_raster(data = relief, aes(x = x, y  = y, alpha = -raytraced), 
                      fill = "gray30",  show.legend = F) +
  ggplot2::geom_sf(data = patch_sjtsk, aes(colour = colour,
                                            fill = colour)) +
  ggplot2::geom_sf(data = krivoklatsko_sjtsk, aes(colour = "black",
                                           fill = NA)) +
  coord_sf(xlim = st_coordinates(bbox_new)[c(1,2),1], # min & max of x values
           ylim = st_coordinates(bbox_new)[c(2,3),2]) + # min & max of y values
  ggplot2::scale_y_continuous(expand = expand_scale(mult = c(0, 0))) +
  ggplot2::scale_x_continuous(expand = expand_scale(mult = c(0, 0))) +
  ggplot2::scale_fill_manual(values = c("black", NA)) +
  ggplot2::scale_colour_manual(values = c("black", "black")) +
  ggplot2::theme_bw() +
  ggplot2::theme(legend.position = "none")
study_map

# Retrieve estimated parameter values ----
a0_m <- sumr[1,2]
a1_m <- sumr[2,2]
a2_m <- sumr[3,2]
b0_m <- sumr[4,2]
b1_m <- sumr[5,2]

# Occupancy simulation ----
no <- nrow(patch_sjtsk)
s <- matrix(0, no, tmax)
s[,1] <- as.numeric(patch_sjtsk$INIT_OCC)
A <- patch_sjtsk %>% dplyr::pull(SHAPEAREA)/10000
P <- patch_sjtsk %>% dplyr::pull(PLANT)

for(t in 1:(tmax-1)) {
  for(i in 1:no) {
    if(s[i,t] == 1) {
      # Persistance
      l <- a0_m+a1_m*P[i]+a2_m*A[i]
      q <- exp(l)/(1+exp(l))
      s[i,t+1] <- rbinom(1,1,q)
    } else {
      # Colonisation
      pj <- rep(0,no)
      for(j in 1:no) {
        l <- b0_m+b1_m*dist[j,i]
        pj[j] <- 1-exp(l)/(1+exp(l))*s[j,t]
      }
      p <- 1-prod(pj)*0.999999
      s[i, t+1] <- rbinom(1,1,p)
    }
  }
}

Oc <- s[,tmax]
si <- s*NA
si[,tmax] <- as.numeric(patch_sjtsk$INIT_OCC)

# Calculate mean occupancy vector 
mean_occ <- rowSums(s)/tmax
