#coding: utf-8

while l = ARGF.gets
  ary = l.chomp.split()

  # skip the comment
  next if ary.first =~ /^#/
  # skip empty line
  next if ary.empty?

  # output and skip the option
  unless ary.first =~ /^\d/
    puts l
    next
  end

  # format and output string
  ary.map! { |e| e.to_f }
  puts "%1.6e %1.6e %1.6e" % ary
end
