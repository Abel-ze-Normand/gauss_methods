$LOAD_PATH << '.'
require 'gauss.rb'
require 'matr_operations.rb'
path = 'matrix.txt'
puts 'Входные данные'
a = Gauss::init path
Gauss::print_matrix a

puts 'Стандартный метод Гаусса'
x = []
x = Gauss::gauss_main_method a

Gauss::print_answer x
err = Gauss::errors path, x

Gauss::norms err

puts '--------------------------------------'
a = Gauss::init path
puts 'Строгий метод Гаусса'
x = Gauss::gauss_strict a

Gauss::print_answer x
err = Gauss::errors path, x

Gauss::norms err

puts '--------------------------------------'
a = Gauss::init path
puts 'LU разложение'
x = Gauss::lu_decomposition a

Gauss::print_answer x
err = Gauss::errors path, x

Gauss::norms err

puts '--------------------------------------'
det = Gauss::det a
puts "Детерминанат А = #{det}"

puts '--------------------------------------'
Gauss::inverse_matrix a


puts '--------------------------------------'

puts "Число обусловленности для различных матричных норм:"
AdvMatr::matrix_norms path

puts '--------------------------------------'
puts "Внесем погрешность 0.01 в правую часть системы"

a = Gauss::init path
av = Gauss::init path
av.each do |arr|
  arr[arr.length - 1] += 0.01
end

x = []
x = Gauss::gauss_main_method a
xv = []
xv = Gauss::gauss_main_method av
puts "Ответ исходной системы и с возмущенной правой частью\n"
Gauss::print_answer x
Gauss::print_answer xv

print "\nОтносительная погрешность:"
norma = x.inject {|memo, item| memo + item * item}
normav = xv.inject {|memo, item| memo +item * item}

norma = Math.sqrt norma
normav = Math.sqrt normav

delta = normav / norma

print "#{format("%.6f", 1 - delta)}\n"
