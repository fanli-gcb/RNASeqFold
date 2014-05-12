#! /usr/bin/ruby

# average pairing probability across every k samples after a burn-in period of MCMC steps
#	struct	ll	u	v	iternum

# NOTE: modified (05/30/2013) to accept per-base digestion probabilities

if ARGV.size == 4
	fn, type, burn_in, k = ARGV
	if !(type == "s" || type == "d")
		puts "type must be either s (for structure) or d (for digestion rates)"
		exit -1
	end
	burn_in = burn_in.to_i
	k = k.to_i
else
	puts "USAGE: #{$0} fn type burn_in k"
	exit -1
end

fp = File.new(fn, "r")
struct, junk = fp.gets.chomp.split(/\t/)
fp.close
n = struct.length

probs = Array.new(n, 0)
rates = Array.new	# [u, v]
1.upto(2) do 
	rates << Array.new(4, 0)
end
count = 0
store_count = 0
i = 0

File.open(fn).each_line do |line|
	count += 1
	line.chomp!
	if count <= burn_in
		next
	end
	
	if i % k == 0
		# store this iteration
		struct, ll, u_str, v_str, r, iternum = line.split(/\t/)
		# structure
		for j in 0..(n-1) do
			if (struct[j].chr == "(" || struct[j].chr == ")")
				probs[j] += 1
			end
		end
		# digestion rates
		u_arr = u_str.chomp.split(",")
		v_arr = v_str.chomp.split(",")
		0.upto(3) do |j|
			rates[0][j] += u_arr[j].to_f
			rates[1][j] += v_arr[j].to_f
		end
		store_count += 1
	end
	i += 1
end

# structure
if type == "s"
	# now divide total count of paired statuses by the number of stored samples
	for j in 0..(n-1) do
		probs[j] = probs[j] * 1.0 / store_count
		puts "#{probs[j]}"
	end
elsif type == "d"
	rates.each { |rate|
		0.upto(3) do |j|
			rate[j] = rate[j] * 1.0 / store_count
		end
		puts rate.join(",")
	}
else
	puts "how did this happen"
	exit -1
end

