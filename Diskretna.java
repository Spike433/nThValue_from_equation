package labos;


import java.util.ArrayList;
import java.util.Scanner;


/******************************************************************************
 *  Compilation:  javac Matrix.java
 *  Execution:    java Matrix
 *
 *  A bare-bones collection of static methods for manipulating
 *  matrices.
 *
 ******************************************************************************/

/*
 * 1. class for matrices  https://introcs.cs.princeton.edu/java/22library/Matrix.java.html  
 * 2. matrix inverse      https://www.sanfoundry.com/java-program-find-inverse-matrix/
 * */
 class Matrix {

  

    // return n-by-n identity matrix I
    public static double[][] identity(int n) {
        double[][] a = new double[n][n];
        for (int i = 0; i < n; i++)
            a[i][i] = 1;
        return a;
    }

    // return x^T y
    public static double dot(double[] x, double[] y) {
        if (x.length != y.length) throw new RuntimeException("Illegal vector dimensions.");
        double sum = 0.0;
        for (int i = 0; i < x.length; i++)
            sum += x[i] * y[i];
        return sum;
    }

    // return B = A^T
    public static double[][] transpose(double[][] a) {
        int m = a.length;
        int n = a[0].length;
        double[][] b = new double[n][m];
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++)
                b[j][i] = a[i][j];
        return b;
    }

    // return c = a + b
    public static double[][] add(double[][] a, double[][] b) {
        int m = a.length;
        int n = a[0].length;
        double[][] c = new double[m][n];
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++)
                c[i][j] = a[i][j] + b[i][j];
        return c;
    }

    // return c = a - b
    public static double[][] subtract(double[][] a, double[][] b) {
        int m = a.length;
        int n = a[0].length;
        double[][] c = new double[m][n];
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++)
                c[i][j] = a[i][j] - b[i][j];
        return c;
    }

    // return c = a * b
    public static double[][] multiply(double[][] a, double[][] b) {
        int m1 = a.length;
        int n1 = a[0].length;
        int m2 = b.length;
        int n2 = b[0].length;
        if (n1 != m2) throw new RuntimeException("Illegal matrix dimensions.");
        double[][] c = new double[m1][n2];
        for (int i = 0; i < m1; i++)
            for (int j = 0; j < n2; j++)
                for (int k = 0; k < n1; k++)
                    c[i][j] += a[i][k] * b[k][j];
        return c;
    }
    
    public static void print(double[][] a) {
    	
    	for(int i=0;i<a.length;i++) {
    		for(int j=0;j<a[i].length;j++) {
    			System.out.print(" "+a[i][j]+" ");
    		}
    		System.out.println();
    	}
    }

    // matrix-vector multiplication (y = A * x)
    public static double[] multiply(double[][] a, double[] x) {
        int m = a.length;
        int n = a[0].length;
        if (x.length != n) throw new RuntimeException("Illegal matrix dimensions.");
        double[] y = new double[m];
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++)
                y[i] += a[i][j] * x[j];
        return y;
    }
    


    // vector-matrix multiplication (y = x^T A)
    public static double[] multiply(double[] x, double[][] a) {
        int m = a.length;
        int n = a[0].length;
        if (x.length != m) throw new RuntimeException("Illegal matrix dimensions.");
        double[] y = new double[n];
        for (int j = 0; j < n; j++)
            for (int i = 0; i < m; i++)
                y[j] += a[i][j] * x[i];
        return y;
    }
  //added from  https://www.sanfoundry.com/java-program-find-inverse-matrix/
    public static double[][] invert(double a[][]) 
    {
        int n = a.length;
        double x[][] = new double[n][n];
        double b[][] = new double[n][n];
        int index[] = new int[n];
        for (int i=0; i<n; ++i) 
            b[i][i] = 1;
 
 // Transform the matrix into an upper triangle
        gaussian(a, index);
 
 // Update the matrix b[i][j] with the ratios stored
        for (int i=0; i<n-1; ++i)
            for (int j=i+1; j<n; ++j)
                for (int k=0; k<n; ++k)
                    b[index[j]][k]
                    	    -= a[index[j]][i]*b[index[i]][k];
 
 // Perform backward substitutions
        for (int i=0; i<n; ++i) 
        {
            x[n-1][i] = b[index[n-1]][i]/a[index[n-1]][n-1];
            for (int j=n-2; j>=0; --j) 
            {
                x[j][i] = b[index[j]][i];
                for (int k=j+1; k<n; ++k) 
                {
                    x[j][i] -= a[index[j]][k]*x[k][i];
                }
                x[j][i] /= a[index[j]][j];
            }
        }
        return x;
    }
    public static void gaussian(double a[][], int index[]) 
    {
        int n = index.length;
        double c[] = new double[n];
 
 // Initialize the index
        for (int i=0; i<n; ++i) 
            index[i] = i;
 
 // Find the rescaling factors, one from each row
        for (int i=0; i<n; ++i) 
        {
            double c1 = 0;
            for (int j=0; j<n; ++j) 
            {
                double c0 = Math.abs(a[i][j]);
                if (c0 > c1) c1 = c0;
            }
            c[i] = c1;
        }
 
 // Search the pivoting element from each column
        int k = 0;
        for (int j=0; j<n-1; ++j) 
        {
            double pi1 = 0;
            for (int i=j; i<n; ++i) 
            {
                double pi0 = Math.abs(a[index[i]][j]);
                pi0 /= c[index[i]];
                if (pi0 > pi1) 
                {
                    pi1 = pi0;
                    k = i;
                }
            }
 
   // Interchange rows according to the pivoting order
            int itmp = index[j];
            index[j] = index[k];
            index[k] = itmp;
            for (int i=j+1; i<n; ++i) 	
            {
                double pj = a[index[i]][j]/a[index[j]][j];
 
 // Record pivoting ratios below the diagonal
                a[index[i]][j] = pj;
 
 // Modify other elements accordingly
                for (int l=j+1; l<n; ++l)
                    a[index[i]][l] -= pj*a[index[j]][l];
            }
        }
    }
}
 

 public class Diskretna {

		public static void main(String[] args) {
			
			
		
			 Scanner input = new Scanner(System.in);
		    	
			 double x0, x1, x2;
			 
			 double a0, a1, a2;

	            System.out.print("Unesite prvo rjesenje x_0 karakteristicne jednadzbe: ");
	            x0 = input.nextDouble();
	            
	            System.out.print("Unesite prvo rjesenje x_1 karakteristicne jednadzbe: ");
	            x1 = input.nextDouble();

	            System.out.print("Unesite prvo rjesenje x_2 karakteristicne jednadzbe: ");
	            x2 = input.nextDouble();

	            System.out.print("Unesite vrijednost nultog clana niza a_0: ");
	            a0 = input.nextDouble();

	            System.out.print("Unesite vrijednost prvog clana niza a_1: ");
	            a1 = input.nextDouble();

	            System.out.print("Unesite vrijednost drugog clana niza a_2: ");
	            a2 = input.nextDouble();

	            System.out.print("Unesite redni broj n trazenog clana niza: ");

	            int x = input.nextInt();
	            
	            double a[][] = {
	            		{Math.pow(x0, 0),Math.pow(x1, 0),Math.pow(x2, 0)},
	            		{Math.pow(x0, 1),Math.pow(x1, 1),Math.pow(x2, 1)},
	            		{Math.pow(x0, 2),Math.pow(x1, 2),Math.pow(x2, 2)}};
	            
	            double c[][]= {
	            		{a0},
	            		
	            		{a1},
	            		
	            		{a2}  // lets make matrix	            	
	            };
	            
	            double b[][]=Matrix.multiply(Matrix.invert(a), c);
	            
	            double rj=b[0][0]*Math.pow(x0, x)+b[1][0]*Math.pow(x1, x)+b[2][0]*Math.pow(x2, x);
	            
	            if((long)rj-rj<1e-6)
	            System.out.println("Vrijednost n-tog clana niza pomocu formule: "+(long )(Math.round(rj * 10d) / 10d));
	            else {
	            	System.out.println("Vrijednost n-tog clana niza pomocu formule: "+(rj * 10d) / 10d);
	            }
		            
				// 2nd way
	            
	            double k0,k1,k2;
	            
	            k0=x0+x1+x2;
	            k1=-(x1*x2)-(x0*x2)-x0*x1;
	            k2=x0*x1*x2;
	            
	            
	           double[] arrInt= new double[x+1];
	           arrInt[0]=a0;
	           arrInt[1]=a1;
	           arrInt[2]=a2;
	            
	            for (int i = 3; i <= x; i++) {
					
	            	arrInt[i]= (k0*arrInt[i-1]+k1*arrInt[i-2]+k2*arrInt[i-3]);
	            	
	            	
				}
	            if((long)arrInt[x]-arrInt[x]<1e-6) {
		            	            
	            System.out.println("Vrijednost n-tog clana niza iz rekurzije: "+(long  )arrInt[x]);
	            
	            }
	            else {
	            	System.out.println("Vrijednost n-tog clana niza iz rekurzije: "+arrInt[x]);
	            }

		}
		
		
	}
 
