public class LU {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		
		double[][] A = {{9,2,1,1},{2,8,-2,1},{-1,-2,7,-2},{1,-1,-2,6}};
		double[] b = {20,16,8,17};

		int A_row = A.length;
		int A_col = A[0].length;
		
		int s = Math.min(A_row, A_col);
		
		/**
		 * L U の形, Aの形により異なる
		 */
		double[][] L;
		double[][] U;
		if(A_row >= A_col)
		{	// Aが縦長の場合 (正方の場合も含める）
			L = new double[A_row][A_col];
			U = new double[A_col][A_col];
		}
		else{
			// Aが横長の場合
			L = new double[A_row][A_row];
			U = new double[A_row][A_col];
		}
		
		double[][] temp_A = copyMatrix(A);
		
		
		for(int i = 0; i < s; i++)
		{
			double[] l_temp = makeTempL(temp_A);
			double[] u_temp = makeTempU(temp_A);
			
			setLMatrix(L, l_temp, i);
			setUMatrix(U, u_temp, i);
			
			double[][] outer = vector_outer(l_temp, u_temp);
			
			double[][] next_A = subMatrix(temp_A, outer);
/*			
			System.out.println(String.format("%d番目のA", i));
			showMatrix(temp_A);
			System.out.println(String.format("%d番目の直積", i));
			showMatrix(outer);
			System.out.println(String.format("%d番目の差", i));
			showMatrix(next_A);
	*/		
			temp_A = reductionCopy(next_A);
			
		}
/*
		System.out.println("L");
		showMatrix(L);
		
		System.out.println("U");
		showMatrix(U);
		System.out.println("LU");
		showMatrix(product(L, U));
		System.out.println("A");
		showMatrix(A);
	*/	
		double[] z = solveLZ_b(L, b);
		double[] x = solveUx_z(U, z);
//		showVector(b);
//		showVector(z);
		showVector(x);
	}
	
	public static void showMatrix(double[][] mat)
	{
		int row = mat.length;
		int col = mat[0].length;
		
		for(int i = 0; i < row; i++)
		{
			System.out.print("|");
			for(int j = 0; j < col; j++)
			{
				System.out.print(String.format("%6.2f", mat[i][j]));
				if(j < col-1)
				{
					System.out.print(", ");
				}
			}
			System.out.println("|");
		}
		
	}
	
	public static void showVector(double[] vec)
	{
		System.out.println(" ");
		int col = vec.length;
		
		
		System.out.print("|");
		for(int i = 0; i < col; i++)
		{
			System.out.print(String.format("%6.2f", vec[i]));
			if(i < col-1)
			{
				System.out.print(", ");
			}
		}
		System.out.println("|");
	}
	
	
	public static double[][] copyMatrix(double[][] inA)
	{
		int row = inA.length;
		int col = inA[0].length;
		
		double[][] temp = new double[row][col];
		
		int t_row = temp.length;
		int t_col = temp[0].length;
		
		for(int i = 0; i < row; i++)
		{
			for(int j = 0; j < col; j++)
			{
				temp[i][j] = inA[i][j];
			}
		}
		return temp;
	}
	
	public static double[][] vector_outer(double[] inA, double[] inB)
	{
		double[][] temp = new double[inA.length][inB.length];
		
		for(int i = 0; i < inA.length; i++)
		{
			for(int j = 0; j < inB.length; j++)
			{
				temp[i][j] = inA[i] * inB[j];
			}
		}
		
		return temp;
	}
	
	public static double[][] subMatrix(double[][] inA, double[][] inB)
	{
		double[][] temp = new double[inA.length][inA[0].length];
		
		for(int i = 0; i < inA.length; i++)
		{
			for(int j = 0; j < inA[0].length; j++)
			{
				temp[i][j] = inA[i][j] - inB[i][j];
			}
		}
		
		return temp;
	}
	
	public static double[] makeTempL(double[][] inA)
	{
		double[] temp_v = new double[inA.length];
		
		for(int i = 0; i < temp_v.length; i++)
		{
			temp_v[i] = inA[i][0] / inA[0][0];
		}
		
		return temp_v;
	}
	
	public static double[] makeTempU(double[][] inA)
	{
		double[] temp_v = new double[inA[0].length];
		
		for(int i = 0; i < temp_v.length; i++)
		{
			temp_v[i] = inA[0][i];
		}
		
		return temp_v;
	}
	
	public static double[][] reductionCopy(double[][] inA)
	{
		int row = inA.length;
		int col = inA[0].length;
		
		double[][] temp = new double[row-1][col-1];
		
		if(row == 1 || col == 1)
		{
			return inA;
		}
		else{
			for(int i = 1; i < row; i++)
			{
				for(int j = 1; j < col; j++)
				{
					temp[i-1][j-1] = inA[i][j];
				}
			}
			
			return temp;
		}
		
	}
	
	public static void setLMatrix(double[][] L, double[] temp_l, int s)
	{
		int row = L.length;
		int l_row = temp_l.length;
		int delta = row - l_row;
		
		for(int i = 0; i < row; i ++)
		{
			if(delta <= i)
			{
				L[i][s] = temp_l[i - delta]; 
			}
			else{
				L[i][s] = 0;
			}
		}
	}
	
	public static void setUMatrix(double[][] U, double[] temp_u, int s)
	{
		int col = U[0].length;
		int u_col = temp_u.length;
		int delta = col - u_col;
		
		for(int j = 0; j < col; j++)
		{
			if(j >= delta)
			{
				U[s][j] = temp_u[j - delta];
			}
			else{
				U[s][j] = 0;
			}
		}
	}
	
	public static double[][] product(double[][] inA, double[][] inB)
	{
		int A_row = inA.length;
		int A_col = inA[0].length;
		int B_row = inB.length;
		int B_col = inB[0].length;
		
		double[][] temp = new double[A_row][B_col];
		
		for(int i = 0; i < A_row; i++)
		{
			for(int j = 0; j < B_col; j++)
			{
				temp[i][j] = 0;
				for(int k = 0; k < B_row; k++)
				{
					temp[i][j] += inA[i][k] * inB[k][j];
				}
			}
		}
		
		return temp;
		
	}
	
	public static double[] solveLZ_b(double[][] L, double[] b)
	{
		int row = L.length;
		int col = L[0].length;
		double[] z = new double[col];
		
//		z[0] = b[0];
//		z[1] = b[1] - L[1][0] * z[0];
//		z[2] = b[2] - L[2][0] * z[0] - L [2][1] * z[1];
		
//		System.out.println(z[0]);
//		System.out.println(z[1]);
//		System.out.println(z[2]);
		for(int i =0; i < col; i++)
		{
			double minus = 0;
			for(int j = 0; j < i; j++)
			{
				minus += L[i][j]*z[j];
			}
			z[i] = b[i] - minus;
		}
		return z;
	}
	
	public static double[] solveUx_z(double[][] U, double[] z)
	{
		int row = U.length;
		int col = U[0].length;
		double[] x = new double[col];
		
//		x[2] = z[2] / U[2][2];
//		x[1] = (z[1] - U[1][2] * x[2]) / U[1][1];
//		x[0] = (z[0] - U[0][2] * x[2] - U[0][1] * x[1]) / U[0][0];
		
		for(int i = col-1 ; i >= 0; i--)
		{
			double minus = 0;
			for(int j = col-1 ; j > i; j--)
			{
				minus += U[i][j] * x[j];
			}
			x[i] = (z[i] - minus) / U[i][i];
		}
		
		return x;
	}


}
