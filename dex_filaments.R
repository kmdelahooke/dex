#DEX_FILAMENTS: convert svg paths into spatstat line segment pattern
#--------------------------------------------------------------------
## arguments: (1) filepath of SVG of the marked up surface as an argument. No large images should be included. 
#(2) retrodeform = TRUE, transform the window based on a retrodeformation factor 
#calculated from the ellipticity of discs, preserving constant area.
## filament id should be left as the default 'pathX'.
## returns a data frame containing the columns x0, x1, y0, y1 
#this can be converted to a spatstat line segment pattern outside the function following:
#df <- dex_filaments(svg, retrodeform = FALSE)
#svg.psp <- psp(df$x0, df$y0, df$x1, df$y1, window = window)


#PACKAGES
if (!require(XML)) install.packages('XML')
if (!require(methods)) install.packages('methods')
if (!require(spatstat)) install.packages('spatstat')
if (!require(spatstat)) install.packages('plyr')
if (!require(spatstat)) install.packages('rlist')
library(XML)
library(methods)
library(spatstat)
library(plyr)
library(rlist)

dex_filaments <- function(filepath, retrodeform = FALSE) {
  #read in SVG
  p <- xmlParse(filepath)
  d <- xmlRoot(p)
  h <- as.numeric(gsub("mm", "", xmlGetAttr(d, "height")))

  #find filaments
  ifelse(length(d[["g"]]["path"]) == 0, {p <- d["path"]}, {p <- d[["g"]]["path"]})

  find_filaments <- function(x) {
    id <- xmlGetAttr(x, "id")
    d <- xmlGetAttr(x, "d")
    if (isFALSE(grepl("z", d)) & grepl("path", id)) {
      return(x)
    }
  }
  p <- list.clean(p, function(x) {length(x) == 0})
  pf <- lapply(p, find_filaments)
  pf <- pf[!sapply(pf, is.null)]

  #extract and split on commands
  fils <- lapply(pf, xmlGetAttr, "d")
  fils <- lapply(fils, strsplit, "(?<=.) (?=[[:alpha:]])", perl = TRUE)
  fils <- lapply(fils, unlist)
  fils <- lapply(fils, strsplit, " ")

  # find number of coords for each command sequence
  length_c <- function(x) {
    cmd <- x[1]
    length <- ifelse(tolower(cmd) == "c", (length(x) - 1) / 3, ifelse(tolower(cmd) == "s", (length(x) - 1) / 2, length(x) - 1))
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
            y1[(l[i - 1] + 1):l[i]] <- cumsum(c(y1[j], x[seq(2, length(x), 2)]))[-1]},
            
      S = { x <- unlist(lapply(x[seq(2, length(x), 2)], strsplit, ","))
            x1[(l[i - 1] + 1):l[i]] <- as.numeric(x[seq(1, length(x), 2)])
            y1[(l[i - 1] + 1):l[i]] <- as.numeric(x[seq(2, length(x), 2)])},

      s = { x <- unlist(lapply(x[seq(2, length(x), 2)], strsplit, ","))
            x1[(l[i - 1] + 1):l[i]] <- cumsum(c(x1[j], x[seq(1, length(x), 2)]))[-1]
            y1[(l[i - 1] + 1):l[i]] <- cumsum(c(y1[j], x[seq(2, length(x), 2)]))[-1]},

      Q = { x <- unlist(lapply(x[seq(2, length(x), 2)], strsplit, ","))
            x1[(l[i - 1] + 1):l[i]] <- as.numeric(x[seq(1, length(x), 2)])
            y1[(l[i - 1] + 1):l[i]] <- as.numeric(x[seq(2, length(x), 2)])},

      q = { x <- unlist(lapply(x[seq(2, length(x), 2)], strsplit, ","))
            x1[(l[i - 1] + 1):l[i]] <- cumsum(c(x1[j], x[seq(1, length(x), 2)]))[-1]
            y1[(l[i - 1] + 1):l[i]] <- cumsum(c(y1[j], x[seq(2, length(x), 2)]))[-1]},

      T = { x <- unlist(lapply(x, strsplit, ","))
            x1[(l[i - 1] + 1):l[i]] <- as.numeric(x[seq(1, length(x), 2)])
            y1[(l[i - 1] + 1):l[i]] <- as.numeric(x[seq(2, length(x), 2)])},

      t = { x <- unlist(lapply(x, strsplit, ","))
            x1[(l[i - 1] + 1):l[i]] <- cumsum(c(x1[j], x[seq(1, length(x), 2)]))[-1]
            y1[(l[i - 1] + 1):l[i]] <- cumsum(c(y1[j], x[seq(2, length(x), 2)]))[-1]})
    }
    df <- data.frame(cbind(x = x1, y = y1))
    return(df)
  }

  #get coords for each path
  s <- NULL
  for (i in seq_len(length(fils))) {
      s[[i]] <- get_coords(fils[[i]])
  }

  print(paste('Number of filaments:', length(s)))

 #get transform
  get_transform <- function(x) {
      t <- xmlGetAttr(x, "transform")
      if (!is.null(t)) {
            t <- gsub("[\\(\\)]", "", regmatches(t, gregexpr("\\(.*?\\)", t))[[1]])
            t <- strsplit(t, ",")
            t <- as.numeric(t[[1]])
            }
      return(t)
      }
  t <- lapply(pf, get_transform)

  #apply path transform
    apply_transform <- function(X, t) {
      x <- X$x
      y <- X$y
      #matrix transform
      ifelse(length(t) == 6, {
          x1 <- t[1] * x + t[3] * y + t[5]
          y <- t[2] * x + t[4] * y + t[6]
          y <- -y + h
      },
      #rotate transform
      {ifelse(length(t) == 1, {
          t <- t * (pi / 180)
          x1 <- cos(t) * x - sin(t) * y
          y <- sin(t) * x + cos(t) * y
          y <- -y + h
      },
      {ifelse(length(t) == 2, {
          x1 <- x + t[1]
          y <- y + t[2]
          y <- -y + h
      }, 
      #no transform
      {
          x1 <- x
          y <- -y + h})
      })}
      )
      return(cbind(x = x1, y = y))
    }

  st <- lapply(mapply(apply_transform, s, t), data.frame)

  get_newxy <- function(X) {
    x0 <- X$x[-length(X$x)]
    x1 <- X$x[-1]
    y0 <- X$y[-length(X$y)]
    y1 <- X$y[-1]
    return(cbind(x0, x1, y0, y1))
  }

  f <- lapply(st, get_newxy)
  f <- ldply(f)
  
  #RETRODEFORMATION
  #find discs
  get_spec <- function(x) {
      id <- xmlGetAttr(x, "id")
      id <- strsplit(id, "_")
      id <- id[[1]]
      if(length(id) == 3) {
            return(x)
        }
  }

  ifelse(length(d[["g"]][c("ellipse", "circle")]) == 0, {ge <- d[c("ellipse", "circle")]}, {ge <- d[["g"]][c("ellipse", "circle")]})
  ge <- lapply(ge, get_spec)
  ge <- ge[!sapply(ge, is.null)]
  
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
  if(length(dsc) == 0) { 
      warning("no discs for retrodeformation")
      return(f)
      }

  #get disc transform
  t <- lapply(dsc, get_transform)
  
  #get disc x, y, rx, ry
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

  dscxy <- ldply(data.frame(mapply(apply_transform, dscxy, t)))[, -1]
  per <- ldply(data.frame(mapply(apply_transform, per, t)))[, -1]
  
  #create list for angle calcs
  coords <- NULL
  for (i in seq_len(length(dscxy[, 1]))) {
      x <- c(dscxy[i, 1], per[i, 1])
      y <- c(dscxy[i, 2], per[i, 2])
      coords[[i]] <- data.frame(cbind(x, y))
  }

  #get_angle
  get_angle <- function(X) {
      r <- atan2(y = (X$y[2] - X$y[1]), x = (X$x[2] - X$x[1])) * 180 / pi - 90
      r[r < 0] <- r + 360
      return(r)
  }
  disc_r <- unlist(lapply(coords, get_angle))
  disc_r <- ifelse(disc_r < 180, 90 - disc_r, 270 - disc_r)
  theta <- - mean(disc_r)
  
  #get ellipicity
  lm <- lm(rx ~ 0 + ry, data = dscxy.df)
  R <- unname(lm$coefficients)
  ifelse(retrodeform == TRUE, {
      print(paste("R squared of  disc rx ~ ry:", signif(summary(lm)$adj.r.squared, 2)))
      print(paste("ellipticity, R:", signif(R, 3)))
      print(paste("theta:", signif(theta, 2)))
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

  f.r <- data.frame(cbind(t(mapply(retro, f$x0, f$y0)), t(mapply(retro, f$x1, f$y1))))
  names(f.r)[1:4] <- c("x0", "y0", "x1", "y1")

  ifelse(retrodeform == FALSE, return(f), return(f.r))
}
