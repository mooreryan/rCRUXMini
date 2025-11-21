checkpoint <- function(checkpoint_file, callback) {
  if (file.exists(checkpoint_file)) {
    NULL
  } else {
    callback()
  }
}
