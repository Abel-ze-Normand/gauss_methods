module AdvMatr
  $LOAD_PATH << '.'
  require 'gauss.rb'
  require 'matrix'

  def self.matrix_square path
    if path.class == String
      temp = Gauss::init path
    else
      temp = path
    end
    matr = Matrix[]
    temp.each_with_index do |arr|
      matr = Matrix.rows(matr.to_a << arr[0..(temp.length - 1)])
    end
    return matr
  end

  def self.matrix_norms path
    #matr_norm 1 = maxj(A)
    norma1 = 0.0
    inv_norma1 = 0.0
    matr = matrix_square path

    puts "--------------Useless output-------------"
    stub = Gauss::init path
    inv_matr = matrix_square Gauss::inverse_matrix stub
    puts "---------End of useless output-----------"

    for j in 0..(matr.column_count - 1)
      sum = 0
      sum_inv = 0
      for i in 0..(matr.column_count - 1)
        sum += matr[i, j].abs
        sum_inv += inv_matr[i, j].abs
      end
      if norma1 < sum
        norma1 = sum
      end
      if inv_norma1 < sum_inv
        inv_norma1 = sum_inv
      end
    end

    puts "||v||1 = #{format("%.6f", norma1 * inv_norma1)}"

    normainf = 0
    inv_normainf = 0

    for i in 0..(matr.column_count - 1)
      sum = 0
      sum_inv = 0
      for j in 0..(matr.column_count - 1)
        sum += matr[i, j].abs
        sum_inv += inv_matr[i, j].abs
      end
      if normainf < sum
        normainf = sum
      end
      if inv_normainf < sum_inv
        inv_normainf = sum_inv
      end
    end

    puts "||v||inf = #{format("%.6f", normainf * inv_normainf)}"

    m = matr * matr.transpose
    m_inv = inv_matr * inv_matr.transpose

    v, d, v_inv = m.eigensystem
    v, d_inv, v_inv = m_inv.eigensystem

    norma2 = Math.sqrt d.each(:diagonal).to_a.max
    inv_norma2 = Math.sqrt d_inv.each(:diagonal).to_a.max

    puts "||v||2 = #{format("%.6f", norma2 * inv_norma2)}"
  end


end
