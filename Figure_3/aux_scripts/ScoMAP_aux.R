#readActivity

readActivity <- function(pathToJpeg){
  img <- readJPEG(pathToJpeg)
  # Obtain the dimension
  imgDm <- dim(img)
  # Assign RGB channels to data frame
  imgRGB <- data.frame(
    x = rep(1:imgDm[2], each = imgDm[1]),
    y = rep(imgDm[1]:1, imgDm[2]),
    R = as.vector(img[,,1]),
    G = as.vector(img[,,2]),
    B = as.vector(img[,,3])
  )
  kClusters <- 2
  kMeans <- kmeans(imgRGB[, c("R", "G", "B")], centers = kClusters)
  clusters <- kMeans$cluster
  signal <- names(which(table(clusters) == min(table(clusters))))
  imgRGB <- imgRGB[which(clusters == signal),]
  image <- sapply(strsplit(pathToJpeg, split = "/"), "[", 3)
  enhancer <- sapply(strsplit(image, split = "__"), "[", 1)
  enhancer <- paste0('GM', enhancer)
  Janelia_line <- rep(enhancer, nrow(imgRGB))
  Image <- rep(image, nrow(imgRGB))
  imgRGB <- cbind(imgRGB, Janelia_line, Image)
  rownames(imgRGB) <- paste0(imgRGB[,1], '_', imgRGB[,2])
  return(imgRGB)
}