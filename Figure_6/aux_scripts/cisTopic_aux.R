# Color palette
distinctColorPalette <-function(k) {
  set.seed(123)
  if(packageVersion("scales") >= '1.1.0'){
    ColorSpace <- t(unique(col2rgb(scales::hue_pal(l=85)(2e3))))
  } else {
    ColorSpace <- t(unique(col2rgb(scales::hue_pal(l=60:100)(2e3))))
  }
  km <- kmeans(ColorSpace, k, iter.max=20)
  colors <- rgb(round(km$centers), maxColorValue=255)
  return(colors)
}
