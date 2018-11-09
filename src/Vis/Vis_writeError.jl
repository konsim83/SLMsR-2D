function writeError2file(err :: Array{Float64,1},
                            t_max :: Float64,
                            file_name :: String, 
                            errType ::String)

  t = range(0, stop=t_max, length=length(err))

  # Write files into #PWD/meshfiles folder
  if ~ispath(pwd() * "/data/error_files")
    mkdir("data/error_files")
  end
  path = pwd() * "/data/error_files/"  
  file_name = path * basename(file_name)

  println("Writing error to   $file_name.*  .......")

  # write points
  open(file_name * ".dat", "w") do f
    firstline = "# Error type:" * errType * "\n"
    write(f, firstline)

    secondline = "# T_max: $(t_max)\n"
    write(f, secondline)

    thirdline = "# Steps: $(length(t_max))\n"
    write(f, thirdline)

    # zip several arrays of different type
    writedlm(f, err, " ")

    close(f)
  end

  return  nothing
end