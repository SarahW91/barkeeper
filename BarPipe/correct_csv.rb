cluster_log = ARGV[0]

new_log = []

count = 0
File.open(cluster_log).each do |line|
	if count == 0
		count += 1
		$columns = line.split("\t").length
		new_log << line
	elsif line.split("\t").length == $columns
		new_log << line
	else
		needed = $columns - line.split("\t").length
		line = line.chomp + "NA" + "\tNA" * needed
		new_log << line
	end
end

File.open(cluster_log, "w+") do |f|
	new_log.each do |line|
		f.puts line
	end
end