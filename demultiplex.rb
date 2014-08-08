## de multiplex glioma data from Peter

require 'zlib'

def main
  barcodeFile = ARGV[0]
  fastq = ARGV[1]

  bar = readBarcodes(barcodeFile)

  
  decode(bar, fastq)

end

def decode(bar, fastq)
  barf = {}
  n = 1
  bar.each_key do |k|
    barf[k] = File.new("cell_#{n}.fastq", 'w')
    n += 1
  end
  
  unknown = File.new("unknown.fastq", 'w')

  if fastq.match(/.gz$/) # gzipped
    fio = Zlib::GzipReader.new(File.open(fastq))
  else
    fio = File.new(fastq, 'r')
  end
    
  chunk = 1000000 # 1 million lines each time
  while !fio.eof?  # not end
    i = 0 

    fh = unknown    
    fio.lines.take(chunk).each do |line|
      i += 1


      if i % 4 == 1  ## first line of a read
        code = line.chomp.split(/\s+/)[-1].split(':')[-1]
        if bar.key?(code)
          fh = barf[code]
          bar[code] += 1
        else
          fh = unknown
        end
      end

      fh.puts line
    end
    
    i = 0
  end

  barf.each do |k,v|
    v.close
  end

  unknown.close
  $stderr.puts bar
end


def readBarcodes(f)
  h = {}
  File.new(f, 'r').each do |line|
    if line=~ /^(\S+)/
      h[$1] = 0
    end
  end
  return h
end



main
