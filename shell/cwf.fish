# CWF fish scripts
# Ismail SEZEN
# sezenismail@gmail.com
# 2021-06-28


# Change directory to defined CWF folder
function goto
  set path (settings.py --path $argv[1])
  if test (count $argv) -eq 2
    set path $path/$argv[2]
  end
  cd $path
end
