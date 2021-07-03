# CWF bash scripts
# Ismail SEZEN
# sezenismail@gmail.com
# 2021-06-28


# Change directory to defined CWF folder
function goto() {
  path=$(settings.py --path $1)
  if [ $# -eq 2 ]
  then
    path="$path/$2"
  fi
  cd $path
}
