#merge all demultiplexed fasta files into one (replacement seqtk fasta)
#for each entry in artificial seqtk fasta add a tag (either trnLF+ITS+rpl16 | matk)

require 'bio'

File.open("/data/data2/lara/Barcoding/2195_intermediate_out/my_seqtk_filtered.fasta", "w+") do |out|
	#collect all files, then iterate through them

	plate_info = {}

	regions = ["trnLF", "rpl16", "ITS", "matK"]
	regions.each do |region|
		files = Dir["/data/data2/lara/Barcoding/2195_intermediate_out/#{region}_parsed/*"]
		files.each do |file_name|
			Bio::FlatFile.open(Bio::FastaFormat, file_name).each do |entry|
				definition = entry.definition.split("|")[1]
				if region != "matK"
					plate_info[definition] = "AAAAAAAAGGGGGGGG"
				else
					plate_info[definition] = "GGGGGGGGAAAAAAAA"
				end
			end
		end
	end

	Bio::FlatFile.open(Bio::FastaFormat, "/data/data2/lara/Barcoding/2195_intermediate_out/seqtk_filtered.fasta").each do |entry|
		out.puts ">" + entry.definition
		if plate_info[entry.definition] != nil
			out.puts plate_info[entry.definition] + entry.seq
		else
			out.puts plate_info.values.sample + entry.seq
		end
	end
end

`rm /data/data2/lara/Barcoding/2195_intermediate_out/seqtk_filtered.fasta`
`mv /data/data2/lara/Barcoding/2195_intermediate_out/my_seqtk_filtered.fasta /data/data2/lara/Barcoding/2195_intermediate_out/seqtk_filtered.fasta`