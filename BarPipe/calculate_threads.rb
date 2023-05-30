idle_cpu = `top -bn1 | grep Cpu`.scan(/\d+,\d\sid/)[0].split(" ")[0].to_f
p `top -bn1 | grep Cpu`
p idle_cpu
free_threads = (idle_cpu/100)*112
num_threads = (free_threads*0.8).to_i
while num_threads%4 != 0
	num_threads -= 1
end
puts num_threads