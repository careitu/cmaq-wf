# CWF bash scripts
# Ismail SEZEN
# sezenismail@gmail.com
# 2021-06-28


# Change directory to defined CWF folder
function goto() {
  path=$(settings.py --path $1)
  cd $path
  echo $PWD
}
