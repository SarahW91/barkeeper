# frozen_string_literal: true

module Exceptions
  class InputError < StandardError; end
  class FinishedPrematurely < StandardError; end
end


class ServerJobManager
  def self.calculate_num_threads
    idle_cpu = `top -bn1 | grep Cpu`.scan(/\d+\.\d\sid/)[0].split(' ')[0].to_f
    free_threads = (idle_cpu / 100) * 112
    num_threads = (free_threads * 0.8).to_i
    num_threads -= 1 while num_threads % 4 != 0
  end
end

class Output
  def expand_tpm

  end
end
