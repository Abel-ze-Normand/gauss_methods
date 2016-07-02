module Gauss
  def self.init path
    myfile = File.new(path)
    matrix = []
    while (!(line = myfile.gets).nil?)
      temp = line.split
      temp.map! { |x| x = x.to_f }
      matrix << temp
    end
    return matrix
  end

  def self.print_matrix a
    a.each do |arr|
      arr.each_with_index do |x, index|
        print "#{format("%.6f", x)} "
        if index == arr.length - 2
          print ' | '
        end
      end
      puts
    end
    puts
  end

  def self.gauss_main_method a
    n = a.length - 1

    for i in 0..n
      (n+1).downto(i) do |j|
        a[i][j] = a[i][j] / a[i][i]
      end

      for k in (i+1)..n
        (n+1).downto(i) do |j|
          a[k][j] = a[k][j] - (a[k][i] * a[i][j])
        end
      end
    end

    puts 'Матрица после прямого хода >>'
    print_matrix a

    x = Array.new(n + 1,0)
    x[n] = a[n][n+1]

    (n-1).downto(0) do |i|
      k = n+1
      for j in i..(n-1)
        a[i][k] = a[i][k] - a[i][j+1] * x[j+1]
        x[i] = a[i][k]
      end
    end
    return x
  end

  def self.gauss_strict a
  #TODO Выбор элемента должен происходить более строго:
  #в таком случае, будет отсутствовать деление на малые числа, следовательно,
  #вычислительная стабильность и точность должны возрасти
  n = a.length - 1

  puts 'Входные данные '
  print_matrix a

  refined = Array.new(n+1, false)
  for i in 0..n
    maxindex = 0
    for j in 0..n #Находим индекс наибольшего элемента
      if a[i][j].abs > a[i][maxindex].abs
        maxindex = j
      end
    end

    max = a[i][maxindex]
    a[i].map! { |j| j / max }
    a[i][maxindex] = 1.0

    for j in 0..n
      if i != j
        coef = a[j][maxindex]
        for k in 0..(n+1)
          unless (refined[k])
            a[j][k] -=coef * a[i][k]
          end
        end
      end
    end
    refined[maxindex] = true;
  end

  puts 'Полученная матрица >>'
  print_matrix a

  x = Array.new(n + 1, 0)
  a.each do |subarr|
    index1 = subarr.index(1.0)
    x[index1] = subarr[n+1]
  end

  return x
  end

  def self.print_answer x
    puts
    print "Ответ: "
    x.each do |num|
      print "#{num.round(6)} "
    end
  end

  def self.lu_decomposition target
    #initialization
    n = target.length - 1
    a = []
    target.each do |arr|
      a << arr[0..n]
    end

    b = []
    for i in 0..n
      b << target[i][n+1]
    end

    l = []
    for i in 0..n
      l << Array.new(n+1, 0)
    end

    u = []
    (n+1).times do
      u << Array.new(n+1, 0)
    end
    #end init
    #main module
    for i in 0..n
      for j in 0..n
        if i <= j
          s = a[i][j]
          for k in 0..(i-1)
            s = s - l[i][k]*u[k][j]
          end
          u[i][j] = s
        end

        if i > j
          s = a[i][j]
          for k in 0..(j-1)
            s = s - l[i][k] * u[k][j]
          end
          l[i][j] = s / u[j][j]
        end
      end
    end
    #end decomposition

    puts 'MATRIX L'
    l.each do |arr|
      arr.each do |c|
        print "#{format("%.6f", c)} "
      end
      puts
    end

    puts 'MATRIX U'
    u.each do |arr|
      arr.each do |c|
        print "#{format("%.6f", c)} "
      end
      puts
    end
    puts
    #схема Холецкого

    res = Array.new(n+1, 0)
    for i in 0..n
      res[i] = b[i]
      for k in 0..(i-1)
        res[i] = res[i] - l[i][k] * res[k]
      end
    end

    x = []

    n.downto(0) do |i|
      for k in (i + 1)..n
        res[i] = res[i] - u[i][k] * res[k]
      end
      res[i] = res[i] / u[i][i]
    end

    x = res
    return x
  end

  def self.det matrix
    a = []
    n = matrix.length - 1
    matrix.each do |arr|
      a << arr[0..n]
    end


    for i in 0..n
      for k in (i+1)..n
        (n).downto(i) do |j|
          a[k][j] = a[k][j] - (a[k][i] * a[i][j]) / a[i][i]
        end
      end
    end

    det = 1.0
    for i in 0..n
      det *= a[i][i]
    end

    return det
  end

  def self.inverse_matrix matrix
    a = []
    n = matrix.length - 1
    matrix.each do |arr|
      a << arr[0..n]
    end

    a_inv = []
    (n+1).times do
      a_inv << Array.new(n+1)
    end

    for i in 0..n
      #construct matrix

      temp = []
      for j in 0..n
        temp << a[j][0..n]
        i == j ? temp[j] << 1.0 : temp[j] << 0.0
      end

      #TODO gauss_main method
      #

      xi = gauss_main_method temp
      xi.each_with_index do |item, index|
        a_inv[index][i] = item
      end

    end


    puts 'Обратная матрица'
    a_inv.each do |arr|
      arr.each do |c|
        print "#{format("%.6f", c)} "
      end
      puts
    end
    puts

    return a_inv
  end

  def self.errors path, x
    puts "\nНевязки\n"

    temp = init path

    errors = []
    temp.each_with_index do |eq, i|
      summ = 0
      eq.each_with_index do |a, index|
        if (!x[index].nil?)
          summ += a * x[index]
        end
      end
      ri = temp[i][temp.length] - summ
      puts "r#{i} = #{ri}"
      errors << ri
    end
    return errors
  end

  def self.norms x
    puts "\n\nНормы\n"

    norma1 = 0
    x.each do |num|
      norma1 += num.abs
    end

    puts "||.||1 = #{norma1}"

    normainf = x.max_by { |a| a.abs}

    puts "||.||inf = #{normainf.abs}"

    norma2 = 0
    x.each do |num|
      norma2 += num * num
    end
    norma2 = Math::sqrt norma2

    puts "||.||2 = #{norma2}"
  end

end
