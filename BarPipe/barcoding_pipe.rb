begin	
	##############
	#USER ARGUMENT
	##############

	require 'optparse'
	require 'ostruct'
	require 'pp'

	class OptparseExample

		def self.parse(args)
			options = OpenStruct.new
			options.map = []
			options.taxon = "Embryophyta"
			options.identity = 0.97
			options.quality = 25
			options.primer_mismatches = nil
			options.tag_mismatches = 2
			options.tpm_splitter = nil
			options.webdav = nil
			options.fa = nil
			options.title = nil

			opt_parser = OptionParser.new do |opts|
				opts.banner = "Usage: barcoding_pipe.rb [options]"
				opts.on("-m", "--require TAGPRIMERMAP",
					"Please provide a tagPrimerMap when executing the script") do |map|
					options.map << map
	     		end
				opts.on("-o", "--require OUTPUT",
					"Please provide an OUTPUT directory when executing the script") do |output_directory|
					options.output = output_directory
				end
				opts.on("-d", "--require ID",
					"Please provide an ID.") do |id|
					options.id = id
				end
				opts.on("-t", "--taxon [TAXA]", "Set taxon") do |taxon|
					taxa_translator = {"Embryophyta" => "Embryophyta", "Marchantiophytina" => "Marchantiophyta", "Anthocerotophytina" => "Anthocerotophyta", "Lycopodiophytina" => "Lycopodiopsida", 
						"Equisetophytina" => "Equisetidae", "Filicophytina" => "Polypodiopsida", "Magnoliopsida" => "Magnoliidae"}
					options.taxon = taxa_translator[taxon]
				end
				opts.on("-i" "--identity THRESHOLD", Float, "Set identity THRESHOLD") do |threshold|
					options.identity = threshold
				end
				opts.on("-q" "--quality THRESHOLD", OptionParser::DecimalInteger, "Set quality THRESHOLD") do |threshold|
					options.quality = threshold
				end
				opts.on("-p" "--primer MISMATCHES", OptionParser::DecimalInteger,  "Set primer MISMATCHES") do |p_mis|
					options.primer_mismatches = p_mis
				end
				opts.on("-b" "--barcode_mismatches", OptionParser::DecimalInteger, "Set barcode mismatches") do |t_mis|
					options.tag_mismatches = t_mis
				end
				opts.on("-s" "--tpm SPLITTER", "Provide SPLITTER file for assigning tagPrimerMap") do |splitter|
					options.tpm_splitter = splitter
				end
				opts.on("-w", "--webdav ADDRESS", "Provide the webdav ADDRESS of a tar file containing bam files") do |address|
					options.webdav = address
				end
				opts.on("-f", "--fa FA", "Please provide either a FASTQ/FASTA file or a webdav address when executing the script") do |fa|
					options.fa = fa
	     		end
	     		opts.on("-j", "--job TITLE", "Please provide a job TITLE when executing the script") do |title|
					options.title = title
	     		end
			end
		    opt_parser.parse!(args)
		    options
		end
	end
	$options = OptparseExample.parse(ARGV)
	$multiple_plates = true if $options["map"].length != 1

	if $options["map"].length > 1 && $options["tpm_splitter"] == nil
		puts "Please provide tags for assigning the according tagPrimerMap."
		abort
	elsif $options["map"].length == 1 && $options["tpm_splitter"] != nil
		puts "You did not provide multiple tagPrimerMaps. Is that correct? [Y|N]"
		yn = gets.chomp
		if yn =~ /[Nn]/
			puts "Please provide multiple TPMs and retry. Thanks."
			abort
		else
			$multiple_plates = false
		end
	end
	if $options["fa"] == nil && $options["webdav"] == nil
		puts "You did not provide a fastq/fasta file or a webdav address. Please do so when running the script again."
		abort
	elsif $options["fa"] != nil && $options["webdav"] != nil
		puts "Please provide EITHER a webdav address or a fastq/fasta file."
		abort
	end

	########
	#METHODS
	########

	require 'tre-ruby'
	require 'fuzzystringmatch'
	require 'ruby-progressbar'
	require 'parallel'
	require 'csv'

	class String
		include TRE
	end

	class QiimeLike
		@jarow = FuzzyStringMatch::JaroWinkler.create(:native)

		def self.splitLibraries(sequence, seq_definition, rc = "noRC")
			@motif_positions = {}
			@sequence = sequence.upcase
			hits = []
			@barcodes = {}
			@primers = {}
			@seq_definition = seq_definition
			@rc_info = rc

			if qualityFiltering(@sequence)
				if $multiple_plates && $options["webdav"] == nil
					@this_plate = findPlate()
					@this_tpm = $tpm_info[@this_plate]
				elsif $multiple_plates
					@this_plate = $plate
					@this_tpm = $tpm_info[@this_plate]
				else
					@this_plate = $tpm_info.keys[0]
					@this_tpm = $tpm_info[@this_plate]
				end

				barcodeRangeGetter
				if @barcodes.length == 0
					$no_barcode == true
				else
					matchingPrimerPairs
					if @primers.length != 0 || @one_primer_only_pairs.length != 0
						sequences = splitMultiples
						sequences_w_motifs = motifAssignmnent(sequences)
						sequences_w_motifs.each do |seq, motifs|
							if qualityFiltering(seq)
								hits << resultPlease(seq, motifs)
							end
						end

						if @regions.length != @regions.uniq.length #delete all but first occurrence of region
							collect_regions = []
							hits.length.times do |i|
								if !collect_regions.include?(hits[i][0])
									collect_regions << hits[i][0]
								else
									hits.delete_at(i)
								end
							end
						elsif @regions.length > 1
							$more_than_one_hit += 1
						end
					end

					if @barcodes.values.uniq.length > 1
						$multiple_barcodes["#{@seq_definition}|#{@rc_info}"] = @barcodes.values
					end
				end
			end
			hits = [[nil, nil]] if hits.empty?

			hits.reject! {|hit| !qualityFiltering(hit[1])}

			return hits
		end

		def self.findPlate
			$splitter_tags.each do |tpm, tag|
				if (@sequence.afind tag, 0, TRE.fuzziness(0)) != nil
					return tpm
				end
			end
		end

		def self.qualityFiltering(seq_to_check)
			if seq_to_check != nil
				if seq_to_check.length.between?(150,7500)
					if seq_to_check.scan(/(\w)\1{31,}/).length != 0
						$homopolymers += 1
						return false
					end
				else
					$too_long_short += 1
					return false
				end
				return true
			else
				return false
			end
		end

		def self.barcodeRangeGetter
			barcode_ranges = {}
			barcode_array = @this_tpm["barcodes"]
			barcode_array.each do |barcode|
				(@sequence.ascan_r barcode, TRE.fuzziness($options["tag_mismatches"])).each do |range|
					barcode_ranges[range.end] = barcode
					@motif_positions[range.end] = "barcode"
				end
			end
			barcode_ranges.values.uniq.each do |b|
				$barcode_counter[b] += 1
				#$tpm_info[@this_plate]["barcode_counter"][b] += 1
			end

			@barcodes = barcode_ranges.sort_by{|k,v| k}.to_h
		end

		def self.matchingPrimerPairs()
			@primers = {}
			@one_primer_only_pairs = []
			primer_hash = @this_tpm["primers"]
			mismatch_hash = @this_tpm["mismatches"]
			primer_hash.each do |linker_primer, reverse_primers|
				l, r, linker_position = false, false, nil
				(@sequence.ascan_r linker_primer, TRE.fuzziness(mismatch_hash[linker_primer])).each do |l_position|
					linker_position = l_position
					l = true
				end
				if l == true
					reverse_primers.each do |reverse_primer|
						(@sequence.ascan_r reverse_primer, TRE.fuzziness(mismatch_hash[reverse_primer])).each do |reverse_position|
							@motif_positions[linker_position.end] = "linker_primer"
							@primers[linker_position.end] = linker_primer
							@motif_positions[reverse_position.begin] = "reverse_primer"
							@primers[reverse_position.begin] = reverse_primer
							r = true
						end
					end
					if r == false
						@one_primer_only_pairs << [linker_primer, reverse_primers]
					end
				end
			end

			@motif_positions = @motif_positions.sort_by{|k,v| k}.to_h
			@primers = @primers.sort_by{|k,v| k}.to_h

			#filtering
			linker = @motif_positions.values.each_index.select { |i| @motif_positions.values[i] == "linker_primer"}
			reverse = @motif_positions.values.each_index.select { |i| @motif_positions.values[i] == "reverse_primer"}
			if !linker.empty? && !reverse.empty?
				#is there a linker before the first reverse?
				if linker[0] > reverse[0]
					@motif_positions.delete(@motif_positions.keys[reverse[0]])
				end
				#is there a reverse after the last linker?
				if reverse[-1] < linker[-1]
					@motif_positions.delete(@motif_positions.keys[linker[-1]])
				end
			end
		end

		def self.splitMultiples()
			@motif_positions = @motif_positions.sort_by{|k, v| -k}.to_h
			last_motif = [-100, nil]
			splitted_seqs = {}
			@motif_positions.each do |position, motif|
				if motif == "barcode"
					motif_length = @barcodes[position].length
				end

				if motif == "barcode" && motif == last_motif[1] && last_motif[0]-(position+motif_length) >= 0 #no splitting where overlap
					splitted_seqs[(position+1...@sequence.length)] = @sequence[position+1..-1]
					@sequence = @sequence[0..position]
				end
				last_motif = [position, motif]
			end

			if splitted_seqs.length == 0
				splitted_seqs = {(0...@sequence.length) => @sequence}
			end

			return splitted_seqs
		end

		def self.motifAssignmnent(seqs)
			sequences = {}
			seqs.each do |pos, seq|
				sequences[seq] = {}
				fitting_bars = []
				ambiguous = false
				@barcodes.each do |b_pos, barcode|
					if b_pos.between?(pos.begin, pos.end+1)
						fitting_bars << barcode
					end
				end
				if fitting_bars.length == 1
					sequences[seq]["barcode"] = fitting_bars[0]
				elsif fitting_bars.length > 1 #more than one barcode found
					bar_counts = Hash.new 0
					fitting_bars.each do |bar|
					  bar_counts[bar] += 1
					end

					if bar_counts.values.count(2) == 1 || fitting_bars.uniq.length == 1
						bar_counts.each { |k,v| sequences[seq]["barcode"] = k if v == 2 }
					else #more than one unique barcode found
						all_motifs = @motif_positions.values
						all_pos = @motif_positions.keys
						all_motifs.length.times do |i|
							if all_pos[i].between?(pos.begin, pos.end+1) && all_motifs[i+3] != nil
								if all_motifs[i] == "barcode" && all_motifs[i+1] == "reverse_primer" && all_motifs[i+2] == "linker_primer" && all_motifs[i+3] == "barcode"
									if @barcodes[all_pos[i]] == @barcodes[all_pos[i+3]]
										sequences[seq]["barcode"] = @barcodes[all_pos[i]]
										break
									end
								end
							end
						end
						if sequences[seq]["barcode"] == nil
							ambiguous = true
						end
					end
				end

				check_for_ambiguous = {"linker_primer" => false, "reverse_primer" => false}
				@primers.each do |p_pos, primer|
					if p_pos.between?(pos.begin, pos.end)
						if @motif_positions[p_pos] == "linker_primer"
							if sequences[seq]["linker_primer"] != nil
								check_for_ambiguous["linker_primer"] = true
							end
							sequences[seq]["linker_primer"] = primer
							sequences[seq]["trunc_start"] = p_pos
							@trunaction_start = p_pos
						elsif @motif_positions[p_pos] == "reverse_primer"
							if sequences[seq]["reverse_primer"] != nil
								check_for_ambiguous["reverse_primer"] = true
							end
							sequences[seq]["reverse_primer"] = primer
							sequences[seq]["trunc_end"] = p_pos
						end
					end
				end

				if (check_for_ambiguous.values.uniq.length == 1 && check_for_ambiguous.values.uniq[0] == true) || ambiguous == true
					File.open("#{$user_output_dir}/ambiguous_seqs.fasta", "a+") do |out|
						out.puts ">#{@seq_definition}|#{@rc_info}|#{pos}|#{seqs.length}|#{@motif_positions.sort_by{|k,v| k}}|#{@primers.sort_by{|k,v| k}}"
						out.puts seq
					end
					$ambiguous_seqs_counter += 1
					sequences[seq] = nil
				end
			end
			return sequences
		end

		def self.resultPlease(sequence, info)
			if info != nil
				barcode = info["barcode"]
				@regions = []
				description_link = @this_tpm["description_link"]
				if info.values.length == 5
					linker_primer = info["linker_primer"]
					reverse_primer = info["reverse_primer"]
					truncation_range = (info["trunc_start"]+1...info["trunc_end"])
					id_description = description_link[[linker_primer, reverse_primer, barcode]]
					if id_description != nil
						@regions << "#{id_description.split("_")[0]}#{id_description.split("_")[-1]}"
						truncated_sequence = sequence[truncation_range]
						return [id_description, truncated_sequence] if truncated_sequence != nil && !truncated_sequence.empty?
					end
				elsif @one_primer_only_pairs.length != 0 #quality filtering bc of one primer only
					$no_second_primer_counter += 1
					@one_primer_only_pairs.each do |linker, reverses|
						reverses.each do |reverse|
							id_description = description_link[[linker, reverse, barcode]]
							$one_primer_only[id_description] = 0 if !$one_primer_only.keys.include?(id_description)
							$one_primer_only[id_description] += 1
							id_description = "no match" if id_description == nil
							File.open("#{$user_output_dir}/incomplete_seqs.fasta", "a+") do |out|
								out.puts ">#{id_description}|#{@seq_definition}|#{@rc_info}"
								out.puts sequence
							end
						end
					end
				end
			end
			return [nil, nil]
		end
	end

	class MyNokogiri
		def self.familyAndGenus(taxonomic_info)
			rgx = /[^sub]family/
			rgx2 = /[^sub]genus/
			family, genus = "NA", "NA"
			rank_info = taxonomic_info.xpath('//LineageEx//Rank').to_a
			scientific_name = taxonomic_info.xpath('//LineageEx//ScientificName').to_a
			(0...rank_info.length).each do |i|
				rank = rank_info[i].to_s
				if rank.match(rgx) != nil
					family = scientific_name[i].content
				elsif rank.include?("genus")
					genus = scientific_name[i].content
					genus = genus.split[0] if genus.include?("sub")
				end
			end
			if family == "NA" && genus == "NA"
				family = scientific_name[scientific_name.length-1].content
			end
			return "\t#{family}\t#{genus}"
		end

		def self.lineage(id)
			taxonomic_info, lineage_info = nil, nil
			while lineage_info == nil do
				taxonomic_info = Nokogiri::XML(Bio::NCBI::REST::EFetch.taxonomy(id, "xml"))
				lineage_info = taxonomic_info.at_xpath('//Lineage')
			end
			lineage = lineage_info.content
			if lineage.include?($options["taxon"])
				lineage_to_print = $options["taxon"]
				addon = familyAndGenus(taxonomic_info)
				if addon == "\tNA\tNA"
					addon = lineage.split(";")[-1].strip
				end
				lineage_to_print = lineage_to_print + addon
			elsif lineage.include?("Fungi")
				lineage_to_print = "Fungi\tNA\tNA"
			elsif lineage.include?("Prokarya")
				lineage_to_print = "Prokarya\tNA\tNA"
			else
				lineage_to_print = "#{lineage}\tNA\tNA"
			end
			if lineage_to_print == nil || lineage_to_print.empty?
				lineage_to_print = "NA\tNA\tNA"
			end
			return lineage_to_print
		end
	end

	def calculateMismatches(primer, allowed_mismatches)
		if allowed_mismatches == nil || allowed_mismatches > 0.75 * primer.length
			wobbles_2_count = primer.scan(/[SWYRKM]/).length
			wobbles_3_count = primer.scan(/[BVDH]/).length
			mismatches = 1 + wobbles_2_count*2 + wobbles_3_count*3
		else
			mismatches = allowed_mismatches
		end
		return mismatches
	end

	def calculateNumThreads
		idle_cpu = `top -bn1 | grep Cpu`.scan(/\d+\.\d\sid/)[0].split(" ")[0].to_f
		free_threads = (idle_cpu/100)*112
		num_threads = (free_threads*0.8).to_i
		while num_threads%4 != 0
			num_threads -= 1
		end

		return 24
	end

	########
	#PRESETS
	########

	require 'bio'

	job_title = $options["title"]

	if Dir.exist?("#{$options["output"]}#{job_title}_intermediate_out")
		puts "It seems like you have run this job before. Do you want to delete the intermediate results before moving on?"
		answer = gets.chomp
		`rm -r #{$options["output"]}#{job_title}_intermediate_out/*parsed` if Dir.exist?("#{$options["output"]}#{job_title}_intermediate_out/*parsed")
		`rm -r '#{$options["output"]}#{job_title}_intermediate_out'` if answer == "yes"
	end
	if !Dir.exist?("#{$options["output"]}#{job_title}_intermediate_out")
		`mkdir '#{$options["output"]}#{job_title}_intermediate_out'`
	end
	$user_output_dir = "#{$options["output"]}#{job_title}_out"
	`rm -r '#{$user_output_dir}'` if Dir.exist?("#{$user_output_dir}")
	`mkdir '#{$user_output_dir}'`

	########
	#COUNTER
	########

	$homopolymers = 0
	$too_long_short = 0
	$ambiguous_seqs_counter = 0
	$no_second_primer_counter = 0
	$barcode_counter = {}
	no_barcode_counter = 0
	$more_than_one_hit = 0

	#######################
	#SEQUEL DATA PROCESSING
	#######################

	if $options["map"].length > 1
		tpms = []
		$options["map"].each do |tpm|
			tpms << tpm.scan(/(?<!\.|\D)[^._\d]+/)[-1]
		end
	else
		tpms = ["A"]
	end

	if $options["webdav"] != nil
		answer = nil
		filtered_fasta = []
		data_folder = "#{$options["output"]}#{job_title}_data"
		webdav_folder = "#{data_folder}/#{job_title}_webdav"

		if File.file?("#{webdav_folder}.tar")
			puts "Sequel data processing was already done. Want to redo it? [Y|N]"
			answer = gets.chomp.downcase
		end

		fasta_folder = "#{data_folder}/fasta"
		if answer == "y" || answer == nil
			ccs_folder = "#{data_folder}/ccs"
			lima_folder = "#{data_folder}/lima"
			
			#out_tar = "#{data_folder}/#{$options["webdav"].split(/[\/\.]/)[-2].gsub(" ", "_")}.tar"
			lima_file = "#{lima_folder}/limaDemux"
			ccs_file = "#{ccs_folder}/#{job_title}.ccs.bam"
			
			`mkdir '#{data_folder}'`
			`mkdir #{ccs_folder} #{lima_folder} #{fasta_folder} #{webdav_folder}`
			
			if !File.file?("#{webdav_folder}.tar")
				`wget -O #{webdav_folder}.tar --user '#{ENV['WEBDAV_USER']}' --password '#{ENV['WEBDAV_PASSWORD']}' '#{$options["webdav"]}'` #fetch .tar from webDav
				`tar -xvf #{webdav_folder}.tar -C #{webdav_folder} --strip-components 1`
			end

			#search for file named "something.subreads.bam" and "something.subreads.bam.pbi" which are the needed starting/input files
			bam = `find #{webdav_folder} -name [^.]*.subreads.bam`.chomp
			pbi = `find #{webdav_folder} -name [^.]*.subreads.bam.pbi`.chomp

			if bam.empty?
				fa = `find #{webdav_folder} -name [^.]*.fastq`.chomp.split("\n")
				type = "fastq"
				if fa.empty?
					fa = `find #{webdav_folder} -name [^.]*.fasta`.chomp.split("\n")
					type = "fasta"
				end
				if fa.length > 1
					`cat #{fa.join(" ")} > #{fasta_folder}/concatenated_subreads.#{type}`
					fa = "#{fasta_folder}/concatenated_subreads.#{type}"
				else
					fa = fa[0]
				end
				$options["fa"] = fa
			else
				# Convert ~subreads.bam to ccs (needs the presence of the pbi-file, as far as I understood)
				puts "Starting CCS"
				if !File.file?(ccs_file)
					`ccs --minLength=100 --num-threads=#{calculateNumThreads} --min-passes=2 --min-rq=0.95 #{bam} #{ccs_file}`
					#`ccs --minLength=100 --num-threads=#{calculateNumThreads()} --min-passes=2 --minReadScore=0.6 --minPredictedAccuracy=0.95 #{bam} #{ccs_file}`
				end

				#following step requires as input a fasta with adapters that define the pools of 96-well plates in addition to the RUN-NUMBER.ccs.bam and a corresponding RUN-NUMBER.ccs.bam.pbi file
				if $multiple_plates
					puts "Starting LIMA"
					if `ls -A #{lima_folder}/`.empty?
						`lima --ccs --split-bam-named -K --num-threads #{calculateNumThreads} #{ccs_file} #{$options["tpm_splitter"]} #{lima_file}.bam` #Demultiplex platepools
					end				
			 
				 	#convert ccs2fasta
					lima_files = `find #{lima_file}*.bam`.split
					lima_files.each do |l_file|
						plates = l_file.scan(/(?<=\.)[^\.]*(?=\.)/)[0]
						if plates.split("-")[0] == plates.split("-")[-1]
							`samtools bam2fq #{lima_file}.#{plates}.bam | seqtk seq -A - > #{fasta_folder}/#{plates.split("-")[0]}.fasta` if !File.file?("#{fasta_folder}/#{job_title}#{plates.split("-")[0]}.fasta")
							filtered_fasta << "#{fasta_folder}/#{plates.split("-")[0]}.fasta"
						end
					end
				else
					`samtools bam2fq #{ccs_file} | seqtk seq -A - > #{fasta_folder}/#{job_title}.fasta` if !File.file?("#{fasta_folder}/#{job_title}.fasta")
					filtered_fasta << "#{fasta_folder}/#{job_title}.fasta"
				end
			end
		else
			filtered_fasta = `find #{fasta_folder}/*.fasta`.split
		end
	end

	############
	#TPM CONTROL
	############

	puts "In need of tpm control?"
	if gets.chomp =~ /[Yy]/
		maps = []
		$options["map"].each do |tpm|
			map_content = File.read(tpm)
			map_content.gsub!(/(\r\r\n|\r\n\r\n|\r\n|\s\#SampleID|ForwardPrimerSequence)/, "\r\r\n" => "\n", "\r\n\r\n" => "\n", "\r\n" => "\n", " #SampleID" => "#SampleID", "ForwardPrimerSequence" => "LinkerPrimerSequence")
			if map_content.include?("\r")
				map_content.gsub!("\r", "\n")
			end

			if map_content.include?("Description")
				delete_last, fields = false, 0
				map_content = map_content.split("\n")
				line_count = 1
				map_content.each do |line|
					if line_count == 1
						fields = line.chomp.split("\t").length
					elsif line_count == 2 && line.chomp.split("\t").length != fields
						delete_last = true
					else
						break
					end
					line_count += 1
				end

				if delete_last == true
					curated_map_content = []
					map_content.each do |line|
						if line.chomp.split("\t").length == fields
							curated_map_content << line.chomp.split("\t")[0...-1].join("\t")
						else
							curated_map_content << line
						end
					end
				end
				map_content = curated_map_content.join("\n")
			end

			map_content.gsub!(/\n\t\t\t\t\t\t/, "")

			altered_tpm = "#{tpm.rpartition(".").first}_altered.txt"
			File.open(altered_tpm, "w+") do |f|
				f.puts map_content
			end

      revised_tpm = `curl -F tpm="@#{altered_tpm}" "https://#{ENV['PROJECT_DOMAIN']}/ngs_runs/revised_tpm" -u "#{ENV['API_USER_NAME']}:#{ENV['API_PASSWORD']}"`

			map_content = revised_tpm.gsub(",", "\t")
			map_content.gsub!("\t\t" "\t")

			File.open(altered_tpm, "w+") do |f|
				f.puts map_content
			end

			maps << altered_tpm
		end
		$options["map"] = maps
	end

	#################
	#QUALITYFILTERING
	#################

	if $options["fa"] != nil && `head -n 4 '#{$options["fa"]}'`.scan(/^@.*\n.*\n\+/).length > 0
		`seqtk seq -a -q #{$options["quality"]} '#{$options["fa"]}' > '#{$options["output"]}#{job_title}_intermediate_out/seqtk_filtered.fasta'`
		filtered_fasta = ["#{$options["output"]}#{job_title}_intermediate_out/seqtk_filtered.fasta"]
	elsif filtered_fasta == nil
		filtered_fasta = [$options["fa"]]
	end

	if $options["webdav"] != nil
		seqs_pre_filtering = "NA"
	end

	seqs_post_filtering_seqtk = 0
	filtered_fasta.each do |fasta|
		seqs_post_filtering_seqtk += `grep -c "^>" '#{fasta}'`.to_i
	end

	###############
	#PLATE SPLITTER
	###############

	if $multiple_plates
		$splitter_tags = {}
		Bio::FlatFile.open(Bio::FastaFormat, $options["tpm_splitter"]).each do |entry|
			$splitter_tags[entry.definition.scan(/[0-9]+[0-9A-Za-z]*/)[-1]] = entry.seq
		end
	end

	########################
	#TAGPRIMERMAP PROCESSING
	########################

	$tpm_info = {}
	regions_array = []
	region_id_definition_link = {}

	$options["map"].each do |tpm|
		plate = tpm.scan(/[0-9]+[0-9A-Za-z]*/)[-1]
		$tpm_info[plate] = {"barcodes" => [], "description_link" => {}, "primers" => {}, "mismatches" => {}, "regions" => []} #"barcode_counter" => {}

		CSV.foreach(tpm, :headers => true, :col_sep => "\t") do |row|
			$tpm_info[plate]["regions"] << row["Region"]
			regions_array << row["Region"]
			barcode, reverse_primer = row["BarcodeSequence"], (Bio::Sequence::NA.new(row["ReversePrimer"]).reverse_complement).upcase
			
			if row["LinkerPrimerSequence"]
				linker_primer = row["LinkerPrimerSequence"]
			elsif row["ForwardPrimerSequence"]
				linker_primer = row["ForwardPrimerSequence"]
			end

			$tpm_info[plate]["barcodes"] << barcode
			$barcode_counter[barcode] = 0
			$tpm_info[plate]["description_link"][[linker_primer, reverse_primer, barcode]] = row["Description"]
			
			if !$tpm_info[plate]["primers"].keys.include?(linker_primer)
				$tpm_info[plate]["primers"][linker_primer] = []
				$tpm_info[plate]["mismatches"][linker_primer] = calculateMismatches(linker_primer, $options["primer_mismatches"])
			end
			if !$tpm_info[plate]["primers"][linker_primer].include?(reverse_primer)
				$tpm_info[plate]["primers"][linker_primer] << reverse_primer
				$tpm_info[plate]["mismatches"][reverse_primer] = calculateMismatches(reverse_primer, $options["primer_mismatches"])
			end

			region_id_definition_link[row["#SampleID"]] = row["Description"]
		end
		$tpm_info[plate]["barcodes"].uniq!
		$tpm_info[plate]["regions"].uniq!
		regions_array.uniq!
	end

	#####################################
	#QIIME REPLACEMENT AKA DEMULTIPLEXING
	#####################################

	progressbar = ProgressBar.create( :format => "%a %b\e[1;33m\u{15E7}\e[0m%i %p%% %t",
		:progress_mark  => ' ',
		:remainder_mark => "\u{FF65}",
		:starting_at    => 0,
		:total => seqs_post_filtering_seqtk)

	no_reverse, run = 0, 0
	no_barcode_forward = false
	$one_primer_only, $seq_to_id, $multiple_barcodes = {}, {}, {}

	filtered_fasta = [filtered_fasta] if !filtered_fasta.kind_of?(Array)

	filtered_fasta.each do |fasta|
		$plate = fasta.scan(/[0-9]+[0-9A-Za-z]*/)[-1]
		Bio::FlatFile.open(Bio::FastaFormat, fasta).each do |entry|
			rc = "noRC"
			$no_barcode = false
			progressbar.increment
			progressbar.format = "%a %b\e[1;33m\u{15E7}\e[0m%i %p%% %t" if run % 10 == 0
			progressbar.format = "%a %b\e[31m\u{15E3}\e[0m%i %p%% %t" if run % 20 == 0

			current_sequence = Bio::Sequence::NA.new(entry.seq)
			$seq_to_id[current_sequence.upcase] = entry.definition
			matched_samples = QiimeLike.splitLibraries(current_sequence, entry.definition)

			no_barcode_forward = $no_barcode

			if matched_samples.uniq == [[nil, nil]]
				rc = "RC"
				$seq_to_id[current_sequence.reverse_complement.upcase] = entry.definition	
				matched_samples = QiimeLike.splitLibraries(current_sequence.reverse_complement, entry.definition, "RC")
				matched_samples[0] == "no reverse" if matched_samples.uniq == [[nil, nil]]
				no_barcode_counter += 1 if $no_barcode == true && no_barcode_forward == true
			end

			if matched_samples.uniq != [[nil, nil]]
				matched_samples.each do |sample|		
					if sample[0] == "no reverse"
						no_reverse += 1
					elsif sample[0] != nil && sample[0] != false
						region = sample[0].split("_")[-1]
						if !Dir.exist?("#{$options["output"]}#{job_title}_intermediate_out/#{region}_parsed")
							`mkdir '#{$options["output"]}#{job_title}_intermediate_out/#{region}_parsed'`
						end

						File.open("#{$options["output"]}#{job_title}_intermediate_out/#{region}_parsed/#{sample[0].split("_")[0]}.fasta", "a+") do |out|
							out.puts ">#{sample[0]}|#{entry.definition}|#{rc}"
							out.puts "#{sample[1].upcase}"
						end
					end
				end
			end
			run += 1
		end
	end

	##################
	#USEARCH AND BLAST
	##################

	require 'nokogiri'

	seqname_to_lineage = {}
	cluster_info = {}
	clusters_per_sample = {}
	cluster_nums = {}
	evalues = {}

	regions_array.each do |region|
		usearch_output_dir = "#{$options["output"]}#{job_title}_intermediate_out/#{region}_clusters"
		`mkdir '#{usearch_output_dir}'` if !Dir.exist?("#{usearch_output_dir}")
		blast_output_dir = "#{$options["output"]}#{job_title}_intermediate_out/#{region}_blast"
		`mkdir '#{blast_output_dir}'` if !Dir.exist?("#{blast_output_dir}")

		region_dir = Dir["#{$options["output"]}#{job_title}_intermediate_out/#{region}_parsed/*"]
		progressbar = ProgressBar.create( :format => "%a %b\e[1;33m\u{15E7}\e[0m%i %p%% %t",
			:progress_mark  => ' ',
			:remainder_mark => "\u{FF65}",
			:starting_at    => 0,
			:total => region_dir.length)
		
		region_dir.each do |parsed_file|
			progressbar.increment
      isolate_id = parsed_file.split(/[\/.]/)[-2]

			if !File.exist?("#{usearch_output_dir}/log_#{isolate_id}")
				`usearch -id #{$options["identity"]} -cluster_fast '#{parsed_file}' -quiet -centroids '#{usearch_output_dir}/#{isolate_id}_centroid' -uc '#{usearch_output_dir}/log_#{isolate_id}' -consout '#{usearch_output_dir}/#{isolate_id}_consensus'`
			end
			if !File.exist?("#{blast_output_dir}/#{isolate_id}")
				puts "blastn -num_threads #{calculateNumThreads} -db nt -max_target_seqs 1 -outfmt '6 qseqid sskingdoms sscinames staxids evalue length pident nident mismatch gaps sacc sseqid stitle' -evalue 1e-10 -query '#{usearch_output_dir}/#{isolate_id}_centroid' -out '#{blast_output_dir}/#{isolate_id}'"
				`blastn -num_threads #{calculateNumThreads} -db nt -max_target_seqs 1 -outfmt '6 qseqid sskingdoms sscinames staxids evalue length pident nident mismatch gaps sacc sseqid stitle' -evalue 1e-10 -query '#{usearch_output_dir}/#{isolate_id}_centroid' -out '#{blast_output_dir}/#{isolate_id}'`
			end

			#how many clusters?
			clusters_per_sample["#{region}#{isolate_id}"] = `grep -c "^C" '#{usearch_output_dir}/log_#{isolate_id}'`.to_i

			#retrieve ids
			seq_name = `awk -F"\t" '{print $1}' '#{blast_output_dir}/#{isolate_id}'`.split("\n")
			taxids = `awk -F"\t" '{print $4}' '#{blast_output_dir}/#{isolate_id}'`.split("\n")
			evalue = `awk -F"\t" '{print $5}' '#{blast_output_dir}/#{isolate_id}'`.split("\n")
			seqname_to_id = seq_name.zip(taxids).to_h

			times = *(0...evalue.length)
			times.each do |i|
				evalues[seq_name[i]] = evalue[i]
			end

			seqname_to_id.each do |seqname, id|
				lineage = MyNokogiri.lineage(id)
				seqname_to_lineage[seqname] = lineage
			end
		end

		`cat '#{usearch_output_dir}/log_'* > '#{usearch_output_dir}/clusterlog_#{region}.txt'`

		#retrieve usearch log file information
		record_type = `awk -F"\t" '{print $1}' '#{usearch_output_dir}/clusterlog_#{region}.txt'`.split("\n")
		cluster_number = `awk -F"\t" '{print $2}' '#{usearch_output_dir}/clusterlog_#{region}.txt'`.split("\n")
		cluster_size = `awk -F"\t" '{print $3}' '#{usearch_output_dir}/clusterlog_#{region}.txt'`.split("\n")
		query_sequence = `awk -F"\t" '{print $9}' '#{usearch_output_dir}/clusterlog_#{region}.txt'`.split("\n")
		target_sequence = `awk -F"\t" '{print $10}' '#{usearch_output_dir}/clusterlog_#{region}.txt'`.split("\n")

		times = *(0...record_type.length)
		times.each do |record_index|
			type = record_type[record_index]
			if type == "C"
				centroid_seq = query_sequence[record_index]
        isolate_id = query_sequence[record_index].split("_")[0]
				cluster_info["#{isolate_id}_#{region}_#{cluster_number[record_index]}"] = [cluster_size[record_index], centroid_seq, seqname_to_lineage[query_sequence[record_index]]]
			end

			cluster_nums[query_sequence[record_index]] = cluster_number[record_index]
			cluster_nums[target_sequence[record_index]] = cluster_number[record_index]
		end
	end

	##################################
	#OUTPUT FILE 1: tagPrimerMap (ext)
	##################################

	$options["map"].each do |map|
		plate = map.scan(/[0-9]+[0-9A-Za-z]*/)[-1]
		this_tpm = $tpm_info[plate]
		CSV.open("#{$user_output_dir}/tagPrimerMap_#{plate}_expanded.txt", "w+", col_sep: "\t") do |out|
			csv = CSV.read(map, :headers => true, :col_sep => "\t")
			out << csv.headers + ["HighQualSeqs", "LinkerPrimerOnly", "Clusters", "TotalSequences"]

			csv.each do |row|
				region = row["Region"]
				description = row["Description"]
				id = row["#SampleID"].scan(/(\D+\d+|\D+\z)/)
        isolate_id = id[0][0]
				id_plus_marker = id.join("_")
				
				if File.exist?("#{$options["output"]}#{job_title}_intermediate_out/#{region}_parsed/#{isolate_id}.fasta")
					high_qual = `grep -c "^>" '#{$options["output"]}#{job_title}_intermediate_out/#{region}_parsed/#{isolate_id}.fasta'`.chomp.to_i
				else
					high_qual = 0
				end
				linker_primer_only = $one_primer_only[description]
				linker_primer_only = 0 if linker_primer_only == nil
				clusters = clusters_per_sample["#{region}#{isolate_id}"]
				out << [id_plus_marker] + row[1..-1] + [high_qual, linker_primer_only, clusters, high_qual + linker_primer_only]
			end
		end
	end

	############################ 
	#OUTPUT FILE 2: cluster info
	############################

	#für jedes cluster drei Spalten: no.seqs, centroid seq, taxFilterResult
	File.open("#{$user_output_dir}/#{job_title}_cluster_log.txt", "w+") do |out|
		out.puts "cluster\tno_seqs\tcentroid_seq\tTaxFilterResult\tfamily\tgenus"
		cluster_info.each do |number, value|
			out.puts "#{number}\t#{value[0]}\t#{value[1]}\t#{value[2]}"
		end
	end

	############################ 
	#OUTPUT FILE 3: barcode info
	############################

	File.open("#{$user_output_dir}/barcode_log.txt", "w+") do |out|
		out.puts "barcode\ttimes_found"
		$barcode_counter.sort_by{ |k, v| -v }.each do |barcode, count|
			out.puts "#{barcode}\t#{count}"
		end
	end

	#################################
	#OUTPUT FILE 3: multiple barcodes
	#################################

	#$multiple_barcodes = $multiple_barcodes.reject {|k,v| v <= 1}
	File.open("#{$user_output_dir}/multiple_barcodes.txt", "w+") do |out|
		out.puts "sequence_id\tno_of_barcodes\tfound_barcodes"
		$multiple_barcodes.sort_by{ |k, v| -v.length }.each do |seq_id, barcodes|
			count = barcodes.length
			str = ""
			barcodes.each do |bar|
				str = "#{str}#{bar},"
			end
			out.puts "#{seq_id}\t#{count}\t#{str[0...-1]}"
		end
	end

	################################# 
	#OUTPUT FILE 5: Multifasta/region
	#################################

	info_for_consensus = {}

	regions_array.each do |region|
		input_fasta_dir = Dir["#{$options["output"]}#{job_title}_intermediate_out/#{region}_parsed/*"]
		File.open("#{$user_output_dir}/#{region}.fasta", "w+") do |out|
			input_fasta_dir.each do |input_fasta_file|
				Bio::FlatFile.open(Bio::FastaFormat, input_fasta_file).each do |entry|
					seqname = entry.definition
					splitted_seqname = seqname.split("|")
          isolate_descr, sequence_id, rc_info = splitted_seqname[0], splitted_seqname[1], splitted_seqname[-1]
					cluster_num = cluster_nums[entry.definition]
          isolate_id = isolate_descr.split("_")[0]
					this_cluster_info = cluster_info["#{isolate_id}_#{region}_#{cluster_num}"]
					
					begin
						lineage, centroid_seqid, cluster_size = this_cluster_info[2], this_cluster_info[1], this_cluster_info[0]
					rescue
						lineage, centroid_seqid, cluster_size = "NA", "NA", "NA"
						File.open("#{$user_output_dir}/#{job_title}_stderr.txt", "a+") do |out|
							out.puts "===\nNo cluster information for the following entry:"
							out.puts entry.definition
							out.puts "#{entry.seq}\n"
						end
					end

					if lineage != nil && lineage != "NA"
						if lineage.include?($options["taxon"])
							lineage_to_print = lineage.split.join(",")
						else
							lineage_to_print = lineage.split("\t")[0]
						end
					else
						lineage_to_print = "NA"
					end

					if isolate_descr == $last_descr
						$num += 1
					else
						$num = 0
					end

					$last_descr = isolate_descr
					isolate_descr = isolate_descr + "_#{$num}"

					#obtain cluster size
					if cluster_num != nil
						if seqname == centroid_seqid
							cluster_num = "#{cluster_num}*centroid"
						end
					end

					evalue = evalues[centroid_seqid]

					definition = ">#{isolate_descr}|#{cluster_num}(#{cluster_size})|#{lineage_to_print}|#{evalue}|#{sequence_id}|#{rc_info}"
					out.puts definition
					out.puts entry.seq

					if definition.include?($options["taxon"])
						File.open("#{$user_output_dir}/#{region}_#{$options["taxon"]}.fasta", "a+") do |e_out|
							e_out.puts definition
							e_out.puts entry.seq
						end
					else
						File.open("#{$user_output_dir}/#{region}_non-#{$options["taxon"]}.fasta", "a+") do |no_e_out|
							no_e_out.puts definition
							no_e_out.puts entry.seq
						end
					end

					#$too_long_short += 1 if definition.include?("NA|NA")
					begin
						info_for_consensus["#{isolate_id}#{region}#{cluster_num.split("*")[0]}"] = [cluster_size, evalue, lineage_to_print]
					rescue
						info_for_consensus["#{isolate_id}#{region}#{cluster_num}"] = [cluster_size, evalue, lineage_to_print]
					end
				end
			end
		end
	end

	##############################
	#OUTPUT FILE 4: consensus seqs
	##############################

	regions_array.each do |region|
		consensus_files = Dir["#{$options["output"]}#{job_title}_intermediate_out/#{region}_clusters/*_consensus"]
		consensus_files.each do |consensus_file|
      isolate_id = consensus_file.scan(/(?<=\/)[^\/]+(?=_)/)[-1]
			Bio::FlatFile.open(Bio::FastaFormat, consensus_file).each do |entry|
				cluster_num = entry.definition.scan(/\d+/)[0]
        isolate_definition = region_id_definition_link["#{isolate_id}#{region}"]
				info = info_for_consensus["#{isolate_id}#{region}#{cluster_num}"]
				cluster_size = info[0]
				taxonomy = info[2]
				evalue = info[1]

				if taxonomy.include?($options["taxon"])
					File.open("#{$user_output_dir}/#{region}_consensus_#{$options["taxon"]}.fasta", "a+") do |e_out|
						e_out.puts ">#{isolate_definition}|#{cluster_num}*consensus(#{cluster_size})|#{taxonomy}|#{evalue}"
						e_out.puts entry.seq
					end
				else
					File.open("#{$user_output_dir}/#{region}_consensus_non-#{$options["taxon"]}.fasta", "a+") do |no_e_out|
						no_e_out.puts ">#{isolate_definition}|#{cluster_num}*consensus(#{cluster_size})|#{taxonomy}|#{evalue}"
						no_e_out.puts entry.seq
					end
				end
			end
		end
	end

	########################## 
	#OUTPUT FILE 3: plate info
	##########################

	seqs_pre_filtering = `wc -l '#{$options["fa"]}'`.to_i / 4 if seqs_pre_filtering != "NA"
	seqs_post_filtering_qiime = 0
	regions_array.each do |region|
		seqs_post_filtering_qiime += `grep -c "^>" '#{$user_output_dir}/#{region}.fasta'`.to_i
	end

	File.open("#{$user_output_dir}/#{job_title}_log.txt", "w+") do |out|
		out.puts "Seqs pre-filtering: #{seqs_pre_filtering}"
		out.puts "Seqs post-seqtk-filtering: #{seqs_post_filtering_seqtk}"
		out.puts "Seqs post-qiimelike-filtering (highQual): #{seqs_post_filtering_qiime}"
		out.puts "Seqs with more than one hit: #{$more_than_one_hit}"
		out.puts "Seqs w/o barcode: #{no_barcode_counter}"
		out.puts "Seqs w/o 2nd primer: #{$no_second_primer_counter}"
		out.puts "Ambiguous seqs: #{$ambiguous_seqs_counter}"
		out.puts "Seqs too short/too long: #{$too_long_short}"
		out.puts "Homopolymers: #{$homopolymers}"
	end

	#zip all output files
	Dir.chdir("#{$options["output"]}"){
	  %x[zip -r '#{job_title}_out.zip' '#{job_title}_out']
	}

	if $options["id"] != nil
    results_path = "#{$options["output"]}#{job_title}_out.zip"
    if ENV['REMOTE_SERVER_PATH']
			`curl -X POST "https://#{ENV['PROJECT_DOMAIN']}/ngs_runs/#{id}/import?results_path=#{$options["output"]}#{job_title}_out.zip" -u "#{ENV['API_USER_NAME']}:#{ENV['API_PASSWORD']}"`
    else
      NgsRun.find(id).import(results_path)
    end
	end
rescue => exception
	File.open("#{$user_output_dir}/#{job_title}_stderr.txt", "a+") do |out|
		out.puts "===\nReason for abort:"
		out.puts exception.backtrace
	end
	raise exception
end
