#DEX: convert svg of marked up fossil surface to a dataframe with coordinate, angle, size and description info for each specimen
#-------------------------------------------------------------------------------------------------------------------------------
#takes the filepath of the svg file as an argument.This svg should contain no images.
#retrodeform = TRUE: transform the window based on a retrodeformation factor 
#calculated from the ellipticity of discs, preserving constant area.
#returns a dataframe with the following columns: specimen id, surface, taxon, x, y, rx, ry, StemL, StemW, StemA, FrondL, FrondW, FrondA, description
#warning messages for: (1) specimen has no disc (hence no coordinates)  (2) length of a path (steml, stemw, frondl, frondw) is zero 
#(3) retrotransformation rx/ry regression is non significant

#PACKAGES
if (!require(XML)) install.packages("XML")
if (!require(methods)) install.packages("methods")
if (!require(plyr)) install.packages("plyr")
library(XML)
library(methods)
library(plyr)

dex <- function(filepath, retrodeform = FALSE) {
  #read in SVG
  p <- xmlParse(filepath)
  d <- xmlRoot(p)
  g <- d[[4]]
  h <- as.numeric(gsub("mm", "", xmlGetAttr(d, "height")))

  #filter out non-specimens
  get_spec <- function(x) {
      id <- xmlGetAttr(x, "id")
      id <- strsplit(id, "_")
      id <- id[[1]]
      if(length(id) == 3) {
          return(x)
      }
  }

  #1.ELLIPSE/CIRCLE: x, y coords, disc radii
  ge <- g[c("ellipse", "circle")]
  ge <- lapply(ge, get_spec)
  ge <- ge[!sapply(ge, is.null)]

  #get x, y, rx, ry
  get_discxy <- function(X) {
      x <- as.numeric(xmlGetAttr(X, "cx"))
      y <- as.numeric(xmlGetAttr(X, "cy"))
      ifelse(xmlName(X) == "ellipse", { 
          rx <- as.numeric(xmlGetAttr(X, "rx"))
          ry <- as.numeric(xmlGetAttr(X, "ry"))
          }, {
          rx <- as.numeric(xmlGetAttr(X, "r"))
          ry <- as.numeric(xmlGetAttr(X, "r"))
          })
      return(cbind(x, y, rx, ry))
  }
  e <- lapply(ge, get_discxy)

  #get transformation
  get_transform <- function(x) {
      t <- xmlGetAttr(x, "transform")
      if (!is.null(t)) {
            t <- gsub("[\\(\\)]", "", regmatches(t, gregexpr("\\(.*?\\)", t))[[1]])
            t <- strsplit(t, ",")
            t <- as.numeric(t[[1]])
      }
      return(t)
  }
  t <- lapply(ge, get_transform)

  #apply transformation
  apply_transform <- function(X, t) {
      x <- X[1]
      y <- X[2]
      ifelse(length(t) == 6, {
          x1 <- t[1] * x + t[3] * y + t[5]
          y <- t[2] * x + t[4] * y + t[6]
          y <- -y + h
      },
      {ifelse(length(t) == 1, {
          t <- t * (pi / 180)
          x1 <- cos(t) * x - sin(t) * y
          y <- sin(t) * x + cos(t) * y
          y <- -y + h
      },
      {
          x1 <- x
          y <- -y + h})
      })
      return(c(x1, y))
  }

  xy <- ldply(data.frame(mapply(apply_transform, e, t)))[, -1]
  names(xy)[1:2] <- c("x", "y")

  #combine transformed x, y, rx, ry, ids, desc into one dataframe
  get_id <- function(x) {
      id <- xmlGetAttr(x, "id")
      id <- strsplit(id, "_")[[1]]
      return(id)
  }

  get_desc <- function(x) {
      desc <- toString.XMLNode(x[1]$desc[1]$text)
      desc[desc == "NULL"] <- 0
      return(desc)
  }

  get_xi <- function(x1, x2) {
    xi <- tolower(paste(x1, x2, collapse = "_"))
    return(xi)
  } 

  rxy <- ldply(e)[, - c(1, 2, 3)]
  id <- ldply(lapply(ge, get_id))[, -1]
  desc <- unname(unlist(lapply(ge, get_desc)))
  xi <- unname(mapply(get_xi, id$V1, id$V2))
  surf <- tolower(id[, 1])
  tx <- tolower(gsub("[[:digit:]]+", "", id[, 2]))

  discs <- data.frame(cbind(xy, rxy, xi, surf, tx, desc))

  #write rows for all discs
  s <- NULL
  for (i in seq_len(nrow(discs))) {
      X <- discs[i, ]
      s[[i]] <- c(X$xi, X$surf, X$tx, X$x, X$y, X$rx, X$ry, 0, 0, 0, 0, 0, 0, X$desc, 0, 0, 0, 0, 0, 0, 0, 0)
  }
  s_discs <- ldply(s)

  #2. PATH: lengths, widths, angles
  gp <- g["path"]
  gp <- lapply(gp, get_spec)
  gp <- gp[!sapply(gp, is.null)]
  p <- lapply(gp, xmlGetAttr, "d")

  #split on commands
  p <- lapply(p, strsplit, "(?<=.) (?=[[:alpha:]])", perl = TRUE)
  p <- lapply(p, unlist)
  p <- lapply(p, strsplit, " ")

  # find number of coords for each command sequence
  length_c <- function(x) {
    cmd <- x[1]
    length <- ifelse(tolower(cmd) == "c", (length(x) - 1) / 3, length(x) - 1)
    return(length)
  }

  #get coords per path
  get_coords <- function(X) {
      l <- cumsum(unlist(lapply(X, length_c))) #cumsum of command sequence lengths
      x1 <- NULL
      y1 <- NULL
      for (i in seq_len(length(l))) {
        cmd <- X[[i]][1] #find command
        x <- X[[i]]
        x <- x[-1]
        j <- l[i - 1] #find index of previous

        switch(cmd,
        M = { x <- unlist(lapply(x, strsplit, ","))
              x1[1:l[i]] <- as.numeric(x[seq(1, length(x), 2)])
              y1[1:l[i]] <- as.numeric(x[seq(2, length(x), 2)])},

        m = { x <- unlist(lapply(x, strsplit, ","))
              x1[1:l[1]] <- cumsum(as.numeric(x[seq(1, length(x), 2)]))
              y1[1:l[1]] <- cumsum(as.numeric(x[seq(2, length(x), 2)]))},

        L = { x <- unlist(lapply(x, strsplit, ","))
              x1[(l[i - 1] + 1):l[i]] <- as.numeric(x[seq(1, length(x), 2)])
              y1[(l[i - 1] + 1):l[i]] <- as.numeric(x[seq(2, length(x), 2)])},

        l = { x <- unlist(lapply(x, strsplit, ","))
              x1[(l[i - 1] + 1):l[i]] <- cumsum(c(x1[j], x[seq(1, length(x), 2)]))[-1]
              y1[(l[i - 1] + 1):l[i]] <- cumsum(c(y1[j], x[seq(2, length(x), 2)]))[-1]},

        V = { x1[(l[i - 1] + 1):l[i]] <- rep(x1[j], length(x))
              y1[(l[i - 1] + 1):l[i]] <- as.numeric(x)},

        v = { ym <- cumsum(c(y1[j], as.numeric(x)))[-1]
              x1[(l[i - 1] + 1):l[i]] <- rep(x1[j], length(ym))
              y1[(l[i - 1] + 1):l[i]] <- ym},

        H = { x1[(l[i - 1] + 1):l[i]] <- as.numeric(x)
              y1[(l[i - 1] + 1):l[i]] <- rep(y1[j], length(x))},

        h = { xm <- cumsum(c(x1[j], as.numeric(x)))[-1]
              x1[(l[i - 1] + 1):l[i]] <- xm
              y1[(l[i - 1] + 1):l[i]] <- rep(y1[j], length(xm))},

        C = { x <- unlist(lapply(x[seq(3, length(x), 3)], strsplit, ","))
              x1[(l[i - 1] + 1):l[i]] <- as.numeric(x[seq(1, length(x), 2)])
              y1[(l[i - 1] + 1):l[i]] <- as.numeric(x[seq(2, length(x), 2)])},

        c = { x <- unlist(lapply(x[seq(3, length(x), 3)], strsplit, ","))
              x1[(l[i - 1] + 1):l[i]] <- cumsum(c(x1[j], x[seq(1, length(x), 2)]))[-1]
              y1[(l[i - 1] + 1):l[i]] <- cumsum(c(y1[j], x[seq(2, length(x), 2)]))[-1]})
      }
     #y1 <- - y1 + h
      df <- data.frame(cbind(x = x1, y = y1))
      return(df)
    }

  #get coords for each path
  coords <- NULL
  for (i in seq_len(length(p))) {
  coords[[i]] <- get_coords(p[[i]])
  }

  coords_2 <- function(X) {
    x0 <- X$x[1]
    x1 <- X$x[2]
    y0 <- X$y[1] 
    y1 <- X$y[2] 
    return(cbind(x0, x1, y0, y1))
  }
  coords2 <- ldply(lapply(coords, coords_2))

  #get id for each path
  id <- lapply(gp, get_id)
  id <- ldply(id)[, -1]

  #find distances
  dist_xy <- function(X) {
      dxy <- sqrt((X$y[2] - X$y[1])^2 + (X$x[2] - X$x[1])^2)
      return(dxy)
  }
  dxy <- unlist(lapply(coords, dist_xy))

  #find angles
  get_angle <- function(X) {
    r <- atan2(y = (X$y[2] - X$y[1]), x = (X$x[2] - X$x[1])) * 180 / pi - 90
    r[r < 0] <- r + 360
    return(r)
  }
  angle <- unlist(lapply(coords, get_angle))

  #combine various descriptors
  desc <- unname(unlist(lapply(gp, get_desc)))
  xi <- unname(mapply(get_xi, id$V1, id$V2))
  surf <- tolower(id[, 1])
  tx <- tolower(gsub("[[:digit:]]+", "", id[, 2]))
  type <- tolower(id[, 3])

  dists <- data.frame(cbind(dxy, type, angle, xi, surf, tx, desc, coords2))

  #write rows based on distance type
  s <- NULL
  for (i in seq_len(nrow(dists))) {
        X <- dists[i, ]
        switch(X$type,
        steml = {
              s[[i]] <- c(X$xi, X$surf, X$tx, 0, 0, 0, 0, X$dxy, 0, X$angle, 0, 0, 0, X$desc, X$x0, X$y0, X$x1, X$y1, 0, 0, 0, 0)},
  
        stemw = {
              s[[i]] <- c(X$xi, X$surf, X$tx, 0, 0, 0, 0, 0, X$dxy, 0, 0, 0, 0, X$desc, 0, 0, 0, 0, 0, 0, 0, 0)},

        frondl = {
              s[[i]] <- c(X$xi, X$surf, X$tx, 0, 0, 0, 0, 0, 0, 0, X$dxy, 0, X$angle, X$desc, 0, 0, 0, 0, X$x0, X$y0, X$x1, X$y1)},

        frondw = {
              s[[i]] <- c(X$xi, X$surf, X$tx, 0, 0, 0, 0, 0, 0, 0, 0, X$dxy, 0, X$desc, 0, 0, 0, 0, 0, 0, 0, 0)}
        )
        
  }

  # combine with ellipses and tidy dataframe 
  s <- ldply(s)
  s <- rbind(s_discs, s)
  names(s)[1:22] <- c("id","surface","taxon", "x","y","rx","ry","StemL", "StemW","StemA", "FrondL", "FrondW","FrondA", "desc", "sx0", "sy0", "sx1", "sy1", "fx0", "fy0", "fx1", "fy1")
  s[is.na(s)] <- 0
  s[, 4:13] <- sapply(s[, 4:13], as.numeric)
  s[, 15:22] <- sapply(s[, 15:22], as.numeric)

  #split on id
  ss <- split(s, s$id)

  #columns reduced to non-zero value
  col_same <- function(x) {
        for (i in seq_len(length(x))) {
              if (x[i] > 0) {
                    x[1] <- x[i]
              }
        }
        return(x[1])
  }

  for (i in seq_len(length(ss))) {
        ss[[i]] <- sapply(ss[[i]], col_same)
  }
  ss <- ldply(ss)[, -1]
  ss[, 4:13] <- sapply(ss[, 4:13], as.numeric)
  ss[, 15:22] <- sapply(ss[, 15:22], as.numeric)

  #remove entries without a disc
  nd <- 0
  for (i in seq_len(length(ss$x))) {
        if(ss$x[i] == 0 & ss$y[i] == 0) {
              warning(paste(ss$id[i], "has no disc"))
              nd <- c(nd, i)
      }
  }
  nd <- nd[-1]

  if (length(nd) > 0) {
      ss <- ss[-nd, ]
  }

  #check stem and frond angles are directed away from disc/pt and remove extra cols
  for (i in seq_len(nrow(ss))) {
    y <- ss$y[i]
    if (ss$StemA[i] != 0) {
        sl1 <- sqrt((- ss$sy1[i] + h - y)^2 + (ss$sx1[i] - ss$x[i])^2)
        sl0 <- sqrt((- ss$sy0[i]  + h - y)^2 + (ss$sx0[i] - ss$x[i])^2)
        ifelse(sl1 > sl0, {
           ifelse(ss$StemA[i] < 180,  {ss$StemA[i] <- ss$StemA[i] + 180}, {ss$StemA[i] <- ss$StemA[i] - 180})    
       }, {ss$StemA[i] <- ss$StemA[i]})
    }
  }
  for (i in seq_len(nrow(ss))) {
    y <- ss$y[i]
    if (ss$FrondA[i] != 0) {
        fl1 <- sqrt((- ss$fy1[i] + h - y)^2 + (ss$fx1[i] - ss$x[i])^2)
        fl0 <- sqrt((- ss$fy0[i] + h - y)^2 + (ss$fx0[i] - ss$x[i])^2)
        ifelse(fl1 > fl0, {
           ifelse(ss$FrondA[i] < 180,  {ss$FrondA[i] <- ss$FrondA[i] + 180}, {ss$FrondA[i] <- ss$FrondA[i] - 180})    
       }, {ss$FrondA[i] <- ss$FrondA[i]})
    }
  }

  ss <- ss[, -seq(15, 22, 1)]

  #RETRODEFORMATION
  #find discs
  get_discs <- function(x) {
      id <- xmlGetAttr(x, "id")
      id <- strsplit(id, "_")[[1]]
      t <- tolower(gsub('[[:digit:]]+', '', id[2]))
      if (t == "disc") {
          return(x)
        }
  }

  dsc <- lapply(ge, get_discs)
  dsc <- dsc[!sapply(dsc, is.null)]
  
  #get disc transform
  t <- lapply(dsc, get_transform)
  
  #get disc x, y, rx, ry
  dscxy <- lapply(dsc, get_discxy)
  dscxy.df <- ldply(dscxy)[, -1]
  
  #perimeter point x1, y1
  get_per <- function (X) {
      x1 <- X[1] + X[3]
      y1 <- X[2]
      return (cbind(x1, y1))
      }  
  per <- lapply(dscxy, get_per)
  
  #apply transform
  dscxy <- ldply(data.frame(mapply(apply_transform, dscxy, t)))[, -1]
  per <- ldply(data.frame(mapply(apply_transform, per, t)))[, -1]
  
  #create list for angle calculations
  coords <- NULL
  for (i in seq_len(length(dscxy[, 1]))) {
      x <- c(dscxy[i, 1], per[i, 1])
      y <- c(dscxy[i, 2], per[i, 2])
      coords[[i]] <- data.frame(cbind(x, y))
    }

  #get angle of disc axis
  get_angle <- function(X) {
      r <- atan2(y = (X$y[2] - X$y[1]), x = (X$x[2] - X$x[1])) * 180 / pi - 90
      r[r < 0] <- r + 360
      return(r)
  }
  disc_r <- unlist(lapply(coords, get_angle))
  disc_r <- ifelse(disc_r < 180, disc_r - 90, disc_r - 270) #difference from x-axis
  theta <- mean(disc_r)
  
  #get ellipicity
  lm <- lm(rx ~ 0 + ry, data = dscxy.df)
  R <- unname(lm$coefficients)
  ifelse(retrodeform == TRUE, {
      print(paste("R squared of  disc rx ~ ry:", signif(summary(lm)$adj.r.squared, 2)))
      print(paste("ellipticity, R:", signif(R, 3)))
      print(paste("theta:", signif(theta, 2)))
      print(paste("p-value:", summary(lm)$coefficients[,4]))
  }, print("no retrodeformation applied"))
  
  if(summary(lm)$coefficients[,4] > 0.05) {
      warnings("disc regression not significant")
  }
  
  #retrotransform
  retro <- function(x, y) {
      t <- theta * pi / 180
      x1 <- cos(t) * x - sin(t) * y
      y <- sin(t) * x + cos(t) * y
      x <- 1 / R * x1 * sqrt(R)
      y <- y * sqrt(R)
      x1 <- cos(-t) * x - sin(-t) * y
      y <- sin(-t) * x + cos(-t) * y
      return (cbind(x = x1, y = y))
  }
  points <- data.frame(t(mapply(retro, ss$x, ss$y)))
  ss2 <- ss
  ss2$x <- points$X1
  ss2$y <- points$X2

ifelse (retrodeform == FALSE, return(ss), return(ss2))
}