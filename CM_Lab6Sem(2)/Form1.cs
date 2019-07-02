using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;

namespace CM_Lab6Sem_2_
{
    public partial class Form1 : Form
    {
        public Form1()
        {
            InitializeComponent();
        }

        public double eps, eps_max1 = 0, eps_max2 = 0, t = 0, xx = 0, yy = 0;
        public int o, s1 = 0, s2 = 0;

       
            
        double f(double x, double y, int o)
        {
            if (o == 0)
            {
                if (x == 1 || x == 2 || y == 2 || y == 3)
                {
                    return Math.Sin(Math.PI * x * y);
                }
                else return Math.Pow(Math.PI, 2) * (-Math.Pow(x, 2) * Math.Sin(Math.PI * x * y) - Math.Pow(y, 2) * Math.Sin(Math.PI * x * y));
            }
            if (o == 1)
            {
                if (x == 1)
                    return (y - 2) * (y - 3);
                if (x == 2)
                    return y * (y - 2) * (y - 3);
                if (y == 2)
                    return (x - 1) * (x - 2);
                if (y == 3)
                    return x * (x - 1) * (x - 2);
                else return -Math.Exp(-x * Math.Pow(y, 2));
            }
            else
                return 0;
        }

        double u(double x,double y)
        {
            return Math.Sin(Math.PI * x * y);
        }

        double m1(double x,double y)
        {
            return (y - 2) * (y - 3);
        }

        double m2(double x, double y)
        {
            return y * (y - 2) * (y - 3);
        }

        double m3(double x, double y)
        {
            return (x - 1) * (x - 2);
        }

        double m4(double x, double y)
        {
            return x * (x - 1) * (x - 2);
        }

        double sk(double[,] A,double[,] B,int n,int m)
        {
            double res = 0;
            for(int i = 1; i < n; i++)
                for(int j = 1; j < m; j++)
                {
                    res += A[i, j] * B[i, j];
                }
            return res;
        }

        double skn1(double[,] A, double[,] B, int n1, int m1, int n2, int m2)
        {
            double sum = 0;
            for (int i = 1; i < n1; i++)
                for (int j = 1; j < m1 + m2; j++)
                    sum += A[i, j] * B[i, j];
            for (int i = n1; i < n1 + n2; i++)
                for (int j = 1; j < m1; j++)
                    sum += A[i, j] * B[i, j];
            return sum;
        }

        void Zeidel(int n, int m, double h, double k, double[,] v, double[] x, double[] y, int o, double e, int N)
        {
            int s = 0;
            double eps_max = 0;
            double h2, k2;
            bool w = false;
            double[] r = new double[m];
            double eps_cur = 0, v_old, v_new, a2;
            h2 = (double)1 / (Math.Pow(h, 2));
            k2 = (double)1 / (Math.Pow(k, 2));
            a2 = 2 * (h2 + k2);
            while (!w)
            {
                eps_max = 0;
                for (int j = 1; j < m; j++)
                {
                    for (int i = 1; i < n; i++)
                    {
                        v_old = v[i, j];
                        v_new = h2 * (v[i + 1, j] + v[i - 1, j]) + k2 * (v[i, j + 1] + v[i, j - 1]);
                        v_new = v_new - f(x[i], y[j], o);
                        v_new = (double)v_new / a2;
                        eps_cur = Math.Abs(v_old - v_new);
                        if (eps_cur > eps_max)
                            eps_max = eps_cur;
                        v[i, j] = v_new;
                    }
                   /* for(int i = 0; i < n; i++)
                    {
                        r[j]=v[i,j]*
                    }*/
                }
                s++;
                if (s >= N || eps_max < e) { w = true; }
            }
            if (s1 == 0)
            {
                s1 = s;
                eps_max1 = eps_max;
            }
            else
            {
                s2 = s;
                eps_max2 = eps_max;
            }
        }

        void SGrad(int n, int m, double h, double k, double [,] v, double [,] B, double[] x, double[] y, int o, double e, int N)
        {
            int s = 0;
            double es = 0;
            double h2, k2;
            bool w = false;
            h2 = (double)1 / (Math.Pow(h, 2));
            k2 = (double)1 / (Math.Pow(k, 2));
            double[,] rs = new double[n + 1, m + 1];
            double[,] zs = new double[n + 1, m + 1];
            double[,] fun = new double[n + 1, m + 1];

            for (int i = 1; i < n; i++)
                for (int j = 1; j < m; j++)
                {
                    rs[i, j] = B[i, j];
                    zs[i, j] = rs[i, j];
                }
            while (!w)
            {
                es = 0;
                double sk1 = sk(rs, rs, n, m);

                for (int i = 1; i < n; i++)
                    for (int j = 1; j < m; j++)
                    {
                        double htmp = h2 * (zs[i - 1, j] + zs[i + 1, j]);
                        double ktmp = k2 * (zs[i, j - 1] + zs[i, j + 1]);
                        fun[i, j] = htmp + ktmp - 2 * (h2 + k2) * zs[i, j];
                    }

                double sk2 = sk(fun, zs, n, m);
                double a = sk1 / sk2;

                for (int i = 1; i < n; i++)
                    for (int j = 1; j < m; j++)
                    {
                        double xs = v[i, j];
                        v[i, j] += a * zs[i, j];
                        if (es < Math.Abs(v[i, j] - xs))
                            es = Math.Abs(v[i, j] - xs);
                    }

                for (int i = 1; i < n; i++)
                    for (int j = 1; j < m; j++)
                        rs[i, j] -= a * fun[i, j];

                sk2 = sk(rs, rs, n, m);
                double b = sk2 / sk1;
                if (b < 0.00000001)
                    w = true;

                for (int i = 1; i < n; i++)
                    for (int j = 1; j < m; j++)
                        zs[i, j] = rs[i, j] + b * zs[i, j];
                s++;
                if (es < e || s >= N) { w = true; }
            }
            if (s1 == 0)
            {
                s1 = s;
                eps_max1 = es;
            }
            else
            {
                s2 = s;
                eps_max2 = es;
            }
        }

        void Sgrad2(int n1, int n2, int m1, int m2, double h1, double k1,double h2,double k2, double[,] v, double[,] B, int o, double e, int N)
        {
            int s = 0;
            double es = 0;
            double h12, k12, h22, k22;
            bool w = false;
            h12 = (double)1 / (h1 * h1);
            k12 = (double)1 / (k1 * k1);
            h22 = 1 / (h2 * h2);
            k22 = 1 / (k2 * k2);
            double[,] rs = new double[n1 + n2 + 1, m1 + m2 + 1];
            double[,] zs = new double[n1 + n1 + 1, m1 + m2 + 1];
            double[,] fun = new double[n1 + n2 + 1, m1 + m2 + 1];

            for (int i = 1; i < n1; i++)
                for (int j = 1; j < m1+m2; j++)
                {
                    rs[i, j] = B[i, j];
                    zs[i, j] = rs[i, j];
                }
            ////////////////////////////////
            for (int i = n1; i < n1+n2; i++)
                for (int j = 1; j < m1; j++)
                {
                    rs[i, j] = B[i, j];
                    zs[i, j] = rs[i, j];
                }
            ////////////////////////////////
            while (!w)
            {
                es = 0;
                double sk1 = skn1(rs, rs, n1, m1, n2, m2);

                for (int i = 1; i < n1; i++)
                    for (int j = 1; j < m1; j++)
                    {
                        double htmp = h12 * (zs[i - 1, j] + zs[i + 1, j]);
                        double ktmp = k12 * (zs[i, j - 1] + zs[i, j + 1]);
                        fun[i, j] = htmp + ktmp - 2 * (h12 + k12) * zs[i, j];
                    }
                //сверху
                for(int i = 1; i < n1; i++)
                {
                    double htmp = h12 * (zs[i - 1, m1] + zs[i + 1, m1]);
                    double ktmp = 2 * zs[i, m1 - 1] / (k1 * (k1 + k2)) + 2 * zs[i, m1 + 1] / (k2 * (k1 + k2));
                    fun[i, m1] = htmp + ktmp - 2 * (h12 + 1 / (k1 * k2)) * zs[i, m1];
                }
                //2 зона
                for(int i=1;i<n1;i++)
                    for(int j = m1 + 1; j < m1 + m2; j++)
                    {
                        double htmp = h12 * (zs[i - 1, j] + zs[i + 1, j]);
                        double ktmp = k22 * (zs[i, j - 1] + zs[i, j + 1]);
                        fun[i, j] = htmp + ktmp - 2 * (h12 + k22) * zs[i, j];
                    }
                //вертикаль
                for(int j = 1; j < m1; j++)
                {
                    double htmp = 2 * zs[n1 - 1, j] / (h1 * (h1 + h2)) + 2 * zs[n1 + 1, j] / (h2 * (h1 + h2));
                    double ktmp = k12 * (zs[n1, j - 1] + zs[n1, j + 1]);
                    fun[n1, j] = htmp + ktmp - 2 * (1 / (h1 * h2) + k12) * zs[n1, j];
                }
                for (int i = n1 + 1; i < n1 + n2; i++)
                    for (int j = 1; j < m1; j++)
                    {
                        double htmp = h22 * (zs[i - 1, j] + zs[i + 1, j]);
                        double ktmp = k12 * (zs[i, j - 1] + zs[i, j + 1]);
                        fun[i, j] = htmp + ktmp - 2 * (h22 + k12) * zs[i, j];
                    }

                double sk2 = skn1(fun, zs, n1, m1, n2, m2);
                double a = sk1 / sk2;

                for (int i = 1; i < n1; i++)
                    for (int j = 1; j < m1+m2; j++)
                    {
                        double xs = v[i, j];
                        v[i, j] += a * zs[i, j];
                        if (es < Math.Abs(v[i, j] - xs))
                            es = Math.Abs(v[i, j] - xs);
                    }

                for (int i = n1; i < n1+n2; i++)
                    for (int j = 1; j < m1; j++)
                    {
                        double xs = v[i, j];
                        v[i, j] += a * zs[i, j];
                        if (es < Math.Abs(v[i, j] - xs))
                            es = Math.Abs(v[i, j] - xs);
                    }

                for(int i = 1; i < n1; i++)
                {
                    for(int j = 1; j < m1 + m2; j++)
                    {
                        rs[i, j] -= a * fun[i, j];
                    }
                }

                for(int i = n1; i < n1 + n2; i++)
                {
                    for(int j = 1; j < m1; j++)
                    {
                        rs[i, j] -= a * fun[i, j];
                    }
                }

                sk2 = skn1(rs, rs, n1, m1, n2, m2);

                double b = sk2 / sk1;
                if (b < 0.0000001)
                    w = true;

                for (int i = 1; i < n1; i++)
                    for (int j = 1; j < m1 + m2; j++)
                        zs[i, j] = rs[i, j] + b * zs[i, j];

                for (int i = n1; i < n1+n2; i++)
                    for (int j = 1; j < m1; j++)
                        zs[i, j] = rs[i, j] + b * zs[i, j];
                s++;
                if (es < e || s >= N) { w = true; }
            }
            if (s1 == 0)
            {
                s1 = s;
                eps_max1 = es;
            }
            else
            {
                s2 = s;
                eps_max2 = es;
            }
        }

        private void button1_Click(object sender, EventArgs e)
        {
            o = 1;
            t = 0;
            eps = Convert.ToDouble(textBox3.Text);
            int n = Convert.ToInt32(textBox1.Text);
            int m = Convert.ToInt32(textBox2.Text);
            int p = Convert.ToInt32(textBox4.Text);

            double a = 1, b = 2, c = 2, d = 3;
            double h = (double)1 / n, h2 = (double)1 / (2 * n), h1 = 1 - h, h11 = 1 - h2;
            double k = (double)1 / m, k2 = (double)1 / (2 * m), k1 = 2 - k, k11 = 2 - k2;

            double[,] v = new double[n + 1, m + 1];
            double[,] v1 = new double[2 * n + 1, 2 * m + 1];
            double[,] x0 = new double[n, m];
            double[] x = new double[n + 1];
            double[] y = new double[m + 1];
            double[] x1 = new double[2 * n + 1];
            double[] y1 = new double[2 * m + 1];
            for (int i = 0; i <= m; i++)
            {
                y[i] = c;
                v[0, i] = f(a, c, o);
                v[n, i] = f(b, c, o);
                c += k;
            }
            c = 2;
            for (int i = 0; i <= n; i++)
            {
                x[i] = a;
                v[i, 0] = f(a, c, o);
                v[i, m] = f(a, d, o);
                a += h;
            }
            a = 1;
            for (int i = 1; i < n; i++)
            {
                for (int j = 1; j < m; j++)
                {
                    v[i, j] = 0;
                    x0[i, j] = v[i, j];
                }
            }
            ///////////////////////
            for (int i = 0; i <= 2*m; i++)
            {
                y1[i] = c;
                v1[0, i] = f(a, c, o);
                v1[2*n, i] = f(b, c, o);
                c += k2;
            }
            c = 2;
            for (int i = 0; i <= 2*n; i++)
            {
                x1[i] = a;
                v1[i, 0] = f(a, c, o);
                v1[i, 2*m] = f(a, d, o);
                a += h2;
            }
            a = 1;
            for (int i = 1; i < 2*n; i++)
            {
                for (int j = 1; j < 2*m; j++)
                {
                    v1[i, j] = 0;
                }
            }

            dataGridView1.Columns.Clear();
            dataGridView1.Rows.Clear();

            for (int i = 0; i < n + 2; i++)
            {
                dataGridView1.Columns.Add("0", string.Format((i - 1).ToString()));
                dataGridView1.Columns[0].HeaderCell.Value = "i";
                dataGridView1[i, 0].Value = h1;
                h1 += h;
                for (int j = 0; j < m + 1; j++)
                {
                    if (i == 0)
                    {
                        dataGridView1.Rows.Add();
                    }
                    dataGridView1.Rows[0].HeaderCell.Value = "j";
                    dataGridView1[0, 0].Value = (string)"Y\\X";
                    dataGridView1.Rows[j + 1].HeaderCell.Value = string.Format((j).ToString(), "0");
                }
                dataGridView1[0, i].Value = k1;
                k1 += k;
            }
            
          /*  for (int i = 0; i < 2*n + 2; i++)
            {
                dataGridView1.Columns.Add("0", string.Format((i - 1).ToString()));
                dataGridView1.Columns[0].HeaderCell.Value = "i";
                dataGridView1[i, 0].Value = h11;
                h11 += h2;
                for (int j = 0; j < 2*m + 1; j++)
                {
                    if (i == 0)
                    {
                        dataGridView1.Rows.Add();
                    }
                    dataGridView1.Rows[0].HeaderCell.Value = "j";
                    dataGridView1[0, 0].Value = (string)"Y\\X";
                    dataGridView1.Rows[j + 1].HeaderCell.Value = string.Format((j).ToString(), "0");
                }
                dataGridView1[0, i].Value = k11;
                k11 += k2;
            }*/


            Zeidel(n, m, h, k, v, x, y, o, eps, p);
            Zeidel(2 * n, 2 * m, h2, k2, v1, x1, y1, o, eps, p);

            int i1 = 2, j1 = 2; 
            for (int i = 1; i < n; i++)
            {
                j1 = 2;
                for (int j = 1; j < m; j++)
                {
                    if (t < Math.Abs(v[i, j] - v1[i1,j1]))
                    {
                        t = Math.Abs(v[i, j] - v1[i1, j1]);
                        xx = x[i];
                        yy = y[j];
                    }
                    j1 += 2;
                }
                i1 += 2;
            }
            /*  for (int i = 1; i < n + 2; i++)
              {
                  for (int j = 1; j <= m + 1; j++)
                  {
                      dataGridView1[i, j].Value = v[i - 1, j - 1];
                  }
              }*/
            /*  for (int i = 1; i < 2*n + 2; i++)
              {
                  for (int j = 1; j <= 2*m + 1; j++)
                  {
                      dataGridView1[i, j].Value = v1[i - 1, j - 1];
                  }
              }*/
            int i2 = 2, j2 = 2;
            for (int i = 1; i < n + 2; i++)
            {
                j2 = 2;
                for (int j = 1; j <= m + 1; j++)
                {
                    dataGridView1[i, j].Value = Math.Abs(v1[i2-2, j2-2] - v[i - 1,j - 1]);
                    j2 += 2;
                }
                i2 += 2;
            }
            i2 = 2;j2 = 2;
            chart1.Series[0].Points.Clear();
            chart2.Series[0].Points.Clear();
            chart3.Series[0].Points.Clear();
            chart1.Series[0].Color = Color.Red;
            chart2.Series[0].Color = Color.DarkBlue;
            chart3.Series[0].Color = Color.Black;
            for (int i = 0; i < n + 1; i++)
            {
                for (int j = 0; j < m + 1; j++)
                {
                    chart1.Series[0].Points.AddXY(x[i], v[i, j]);
                    chart2.Series[0].Points.AddXY(x1[i], v1[i, j]);
                    chart3.Series[0].Points.AddXY(x[i], Math.Abs(v[i, j] - v1[2*i,2*j]));
                }
            }
        }

        private void button2_Click(object sender, EventArgs e)
        {
            o = 0;
            t = 0;
            int n = Convert.ToInt32(textBox1.Text);
            int m = Convert.ToInt32(textBox2.Text);
            eps = Convert.ToDouble(textBox3.Text);
            int p = Convert.ToInt32(textBox4.Text);

            double a = 1, b = 2, c = 2, d = 3;
            double h = (double)1 / n, h1 = 1 - h;
            double k = (double)1 / m, k1 = 2 - k;
            double[,] v = new double[n + 1, m + 1];
            double[] x = new double[n + 1];
            double[] y = new double[m + 1];

            for (int i = 0; i <= m; i++)
            {
                y[i] = c;
                v[0, i] = f(a, c, o);
                v[n, i] = f(b, c, o);
                c += k;
            }
            c = 2;
            for (int i = 0; i <= n; i++)
            {
                x[i] = a;
                v[i, 0] = f(a, c, 0);
                v[i, m] = f(a, d, 0);
                a += h;
            }
            a = 1;
            for (int i = 1; i < n; i++)
            {
                for (int j = 1; j < m; j++)
                {
                    v[i, j] = 0;
                }
            }

            dataGridView1.Columns.Clear();
            dataGridView1.Rows.Clear();

            for (int i = 0; i < n + 2; i++)
            {
                dataGridView1.Columns.Add("0", string.Format((i - 1).ToString()));
                dataGridView1.Columns[0].HeaderCell.Value = "i";
                dataGridView1[i, 0].Value = h1;
                h1 += h;
                for (int j = 0; j < m + 1; j++)
                {
                    if (i == 0)
                    {
                        dataGridView1.Rows.Add();
                    }
                    dataGridView1.Rows[0].HeaderCell.Value = "j";
                    dataGridView1[0, 0].Value = (string)"Y\\X";
                    dataGridView1.Rows[j + 1].HeaderCell.Value = string.Format((j).ToString(), "0");
                }
                dataGridView1[0, i].Value = k1;
                k1 += k;
            }
            Zeidel(n, m, h, k, v, x, y, o, eps, p);
            for (int i = 1; i < n + 2; i++)
            {
                for (int j = 1; j <= m + 1; j++)
                {
                    dataGridView1[i, j].Value = Math.Abs(u(x[i - 1], y[j - 1]) - v[i - 1, j - 1]);
                }
            }
            for(int i = 1; i < n; i++)
            {
                for(int j = 1; j < m; j++)
                {
                    if (t < Math.Abs(v[i, j] - u(x[i], y[j])))
                    {
                        t = Math.Abs(v[i, j] - u(x[i], y[j]));
                        xx = x[i];
                        yy = y[j];
                    }
                }
            }
            chart1.Series[0].Points.Clear();
            chart2.Series[0].Points.Clear();
            chart3.Series[0].Points.Clear();
            chart1.Series[0].Color = Color.Red;
            chart2.Series[0].Color = Color.DarkBlue;
            chart3.Series[0].Color = Color.Black;
            for(int i = 0; i < n + 1; i++)
                for(int j = 0; j < m + 1; j++)
                {
                    chart1.Series[0].Points.AddXY(x[i], u(x[i],y[j]));
                    chart2.Series[0].Points.AddXY(x[i], v[i, j]);
                    chart3.Series[0].Points.AddXY(x[i], Math.Abs(v[i, j] - u(x[i], y[j])));
                }
        }

        private void button4_Click(object sender, EventArgs e)
        {
            o = 0;
            t = 0;
            int n = Convert.ToInt32(textBox1.Text);
            int m = Convert.ToInt32(textBox2.Text);
            eps = Convert.ToDouble(textBox3.Text);
            int p = Convert.ToInt32(textBox4.Text);

            double a = 1, b = 2, c = 2, d = 3;
            double h = (double)1 / n, h1 = 1 - h;
            double k = (double)1 / m, k1 = 2 - k;
            double[,] v = new double[n + 1, m + 1];
            double[,] x0 = new double[n + 1, m + 1];
            double[] x = new double[n + 1];
            double[] y = new double[m + 1];

            for (int i = 0; i <= m; i++)
            {
                y[i] = c;
                c += k;
            }
            c = 2;
            for (int i = 0; i <= n; i++)
            {
                x[i] = a;
                a += h;
            }
            a = 1;

            dataGridView1.Columns.Clear();
            dataGridView1.Rows.Clear();

            for (int i = 0; i < n + 2; i++)
            {
                dataGridView1.Columns.Add("0", string.Format((i - 1).ToString()));
                dataGridView1.Columns[0].HeaderCell.Value = "i";
                dataGridView1[i, 0].Value = h1;
                h1 += h;
                for (int j = 0; j < m + 1; j++)
                {
                    if (i == 0)
                    {
                        dataGridView1.Rows.Add();
                    }
                    dataGridView1.Rows[0].HeaderCell.Value = "j";
                    dataGridView1[0, 0].Value = (string)"Y\\X";
                    dataGridView1.Rows[j + 1].HeaderCell.Value = string.Format((j).ToString(), "0");
                }
                dataGridView1[0, i].Value = k1;
                k1 += k;
            }

            for (int i = 0; i < n + 1; i++)
                for (int j = 0; j < m + 1; j++)
                {
                    v[i, j] = 0;
                    x0[i, j] = f(x[i], y[j],o);
                    if (i - 1 == 0)
                        x0[i, j] -= u(x[o], y[j]) / (h * h);
                    if (i + 1 == n)
                        x0[i, j] -= u(x[n], y[j]) / (h * h);
                    if (j - 1 == 0)
                        x0[i, j] -= u(x[i], y[0]) / (k * k);
                    if (j + 1 == m)
                        x0[i, j] -= u(x[i], y[m]) / (k * k);
                }

            SGrad(n, m, h, k, v, x0, x, y, o, eps, p);

            for (int i = 0; i <= m; i++)
            {
                v[0, i] = f(a, c, o);
                v[n, i] = f(b, c, o);
                c += k;
            }
            c = 2;
            for (int i = 0; i <= n; i++)
            {
                v[i, 0] = f(a, c, 0);
                v[i, m] = f(a, d, 0);
                a += h;
            }
            a = 1;

            for (int i = 1; i < n + 2; i++)
            {
                for (int j = 1; j <= m + 1; j++)
                {
                    dataGridView1[i, j].Value = v[i - 1, j - 1];
                }
            }
            for (int i = 1; i < n; i++)
            {
                for (int j = 1; j < m; j++)
                {
                    if (t < Math.Abs(v[i, j] - u(x[i], y[j])))
                    {
                        t = Math.Abs(v[i, j] - u(x[i], y[j]));
                        xx = x[i];
                        yy = y[j];
                    }
                }
            }
        }

        private void button5_Click(object sender, EventArgs e)
        {
            o = 1;
            t = 0;
            eps = Convert.ToDouble(textBox3.Text);
            int n = Convert.ToInt32(textBox1.Text);
            int m = Convert.ToInt32(textBox2.Text);
            int p = Convert.ToInt32(textBox4.Text);

            double a = 1, b = 2, c = 2, d = 3;
            double h = (double)1 / n, h2 = (double)1 / (2 * n), h1 = 1 - h;
            double k = (double)1 / m, k2 = (double)1 / (2 * m), k1 = 2 - k;

            double[,] v = new double[n + 1, m + 1];
            double[,] v1 = new double[2 * n + 1, 2 * m + 1];
            double[,] x0 = new double[n + 1, m + 1];
            double[,] x00 = new double[2 * n + 1, 2 * m + 1];
            double[] x = new double[n + 1];
            double[] y = new double[m + 1];
            double[] x1 = new double[2 * n + 1];
            double[] y1 = new double[2 * m + 1];

            for (int i = 0; i <= m; i++)
            {
                y[i] = c;
                c += k;
            }
            c = 2;
            for (int i = 0; i <= n; i++)
            {
                x[i] = a;
                a += h;
            }
            a = 1;
            ///////////////////////
            for (int i = 0; i <= 2 * m; i++)
            {
                y1[i] = c;
                c += k2;
            }
            c = 2;
            for (int i = 0; i <= 2 * n; i++)
            {
                x1[i] = a;
                a += h2;
            }
            a = 1;

            dataGridView1.Columns.Clear();
            dataGridView1.Rows.Clear();

            for (int i = 0; i < n + 2; i++)
            {
                dataGridView1.Columns.Add("0", string.Format((i - 1).ToString()));
                dataGridView1.Columns[0].HeaderCell.Value = "i";
                dataGridView1[i, 0].Value = h1;
                h1 += h;
                for (int j = 0; j < m + 1; j++)
                {
                    if (i == 0)
                    {
                        dataGridView1.Rows.Add();
                    }
                    dataGridView1.Rows[0].HeaderCell.Value = "j";
                    dataGridView1[0, 0].Value = (string)"Y\\X";
                    dataGridView1.Rows[j + 1].HeaderCell.Value = string.Format((j).ToString(), "0");
                }
                dataGridView1[0, i].Value = k1;
                k1 += k;
            }

            for (int i = 0; i < n + 1; i++)
                for (int j = 0; j < m + 1; j++)
                {
                    v[i, j] = 0;
                    x0[i, j] = f(x[i], y[j], o);
                    if (i - 1 == 0)
                        x0[i, j] -= m1(x[o], y[j]) / (h * h);
                    if (i + 1 == n)
                        x0[i, j] -= m2(x[n], y[j]) / (h * h);
                    if (j - 1 == 0)
                        x0[i, j] -= m3(x[i], y[0]) / (k * k);
                    if (j + 1 == m)
                        x0[i, j] -= m4(x[i], y[m]) / (k * k);
                }

            for (int i = 0; i < 2*n + 1; i++)
                for (int j = 0; j < 2*m + 1; j++)
                {
                    v1[i, j] = 0;
                    x00[i, j] = f(x1[i], y1[j], o);
                    if (i - 1 == 0)
                        x00[i, j] -= m1(x1[o], y1[j]) / (h2 * h2);
                    if (i + 1 == 2 * n)
                        x00[i, j] -= m2(x1[2 * n], y1[j]) / (h2 * h2);
                    if (j - 1 == 0)
                        x00[i, j] -= m3(x1[i], y1[0]) / (k2 * k2);
                    if (j + 1 == 2 * m)
                        x00[i, j] -= m4(x1[i], y1[2 * m]) / (k2 * k2);
                }

            SGrad(n, m, h, k, v, x0, x, y, o, eps, p);
            SGrad(2 * n, 2 * m, h2, k2, v1, x00, x1, y1, o, eps, p);

            for (int i = 0; i <= m; i++)
            {
                v[0, i] = f(a, c, o);
                v[n, i] = f(b, c, o);
                c += k;
            }
            c = 2;
            for (int i = 0; i <= n; i++)
            {
                v[i, 0] = f(a, c, o);
                v[i, m] = f(a, d, o);
                a += h;
            }
            a = 1;

            int i1 = 2, j1 = 2;
            for (int i = 1; i < n; i++)
            {
                j1 = 2;
                for (int j = 1; j < m; j++)
                {
                    if (t < Math.Abs(v[i, j] - v1[i1, j1]))
                    {
                        t = Math.Abs(v[i, j] - v1[i1, j1]);
                        xx = x[i];
                        yy = y[j];
                    }
                    j1 += 2;
                }
                i1 += 2;
            }
            for (int i = 1; i < n + 2; i++)
            {
                for (int j = 1; j <= m + 1; j++)
                {
                    dataGridView1[i, j].Value = v[i - 1, j - 1];
                }
            }
        }

        private void button6_Click(object sender, EventArgs e)
        {
            o = 0;
            t = 0;
            int n = Convert.ToInt32(textBox1.Text);
            int m = Convert.ToInt32(textBox2.Text);
            eps = Convert.ToDouble(textBox3.Text);
            int p = Convert.ToInt32(textBox4.Text);
            int n1, n2, m1, m2;

            n1 = n / 2;
            n2 = n - n1;
            m1 = m / 2;
            m2 = m - m1;

            double a = 1, b = 2, c = 2, d = 3;
            double h =(double) 1 / n, k =(double) 1 / m;
            double H1 = (double)0.5 / n1, H2 = (double)0.5 / n2, h1 = 1 - h;
            double K1 = (double)0.5 / m1, K2 =(double) 0.5 / m2, k1 = 2 - k;
            double[,] v = new double[n + 1, m + 1];
            double[,] x0 = new double[n + 1, m + 1];
            double[,] U = new double[n + 1, m + 1];
            double[] x = new double[n + 1];
            double[] y = new double[m + 1];

            for (int i = 0; i <= m; i++)
            {
                y[i] = c;
                c += k;
            }
            c = 2;

            for (int i = 0; i <= n; i++)
            {
                x[i] = a;
                a += h;
            }
            a = 1;

            dataGridView1.Columns.Clear();
            dataGridView1.Rows.Clear();

            for (int i = 0; i < n + 2; i++)
            {
                dataGridView1.Columns.Add("0", string.Format((i - 1).ToString()));
                dataGridView1.Columns[0].HeaderCell.Value = "i";
                dataGridView1[i, 0].Value = h1;
                h1 += h;
                for (int j = 0; j < m + 1; j++)
                {
                    if (i == 0)
                    {
                        dataGridView1.Rows.Add();
                    }
                    dataGridView1.Rows[0].HeaderCell.Value = "j";
                    dataGridView1[0, 0].Value = (string)"Y\\X";
                    dataGridView1.Rows[j + 1].HeaderCell.Value = string.Format((j).ToString(), "0");
                }
                dataGridView1[0, i].Value = k1;
                k1 += k;
            }

            for (int i = 0; i < n1; i++)
                for (int j = 0; j < m1; j++)
                {
                    v[i, j] = 0;
                    U[i, j] = u(x[i], y[j]);
                    x0[i, j] = f(x[i], y[j], o);
                    if (i - 1 == 0)
                        x0[i, j] -= u(x[0], y[j]) / (H1 * H1);
                    if (j - 1 == 0)
                        x0[i, j] -= u(x[i], y[0]) / (K1 * K1);
                }

            for (int i = 0; i < n1+1; i++)
                for (int j = m1; j < m1+m2+1; j++)
                {
                    v[i, j] = 0;
                    U[i, j] = u(x[i], y[j]);
                    x0[i, j] = f(x[i], y[j], o);
                    if (i - 1 == 0)
                        x0[i, j] -= u(x[0], y[j]) / (H1 * H1);
                    if (j + 1 == m1+m2)
                        x0[i, j] -= u(x[i], y[m1+m2]) / (K2 * K2);
                    if (i + 1 == n1)
                        x0[i, j] -= u(x[n1], y[j]) / (H1 * H1);
                }

            for (int i = n1; i < n1+n2 + 1; i++)
                for (int j = 0; j < m1 + 1; j++)
                {
                    v[i, j] = 0;
                    U[i, j] = u(x[i], y[j]);
                    x0[i, j] = f(x[i], y[j], o);
                    if (i + 1 == n1+n2)
                        x0[i, j] -= u(x[n1+n2], y[j]) / (H2 * H2);
                    if (j - 1 == 0)
                        x0[i, j] -= u(x[i], y[0]) / (K1 * K1);
                    if (j + 1 == m1)
                        x0[i, j] -= u(x[i], y[m1]) / (K1 * K1);
                }

            Sgrad2(n1, n2, m1, m2, H1, K1, H2, K2, v, x0, o, eps, p);

            for (int i = 0; i <= m; i++)
            {
                v[0, i] = u(a, c);
                c += k;
            }
            c = 2;
            for (int i = 0; i <= n; i++)
            {
                v[i, 0] = u(a, c);
                a += h;
            }
            a = 1;
            for (int i = 0; i <= m1; i++)
            {
                v[n, i] = u(a, c);
                c += k;
            }
            c = 2;
            for (int i = 0; i <= n1; i++)
            {
                v[i, m] = u(a, c);
                a += h;
            }
            a = 1;
            double a1 = 1.5, c1 = 2.5;
            for (int i = m1; i < m1+m2; i++)
            {
                v[n1, i] = u(a1, c1);
                c1 += k;
            }
            c1 = 1.5;
            for (int i = n1+1; i < n1+n2; i++)
            {
                v[i, m1] = u(a1, c1);
                a1 += h;
            }
            a1 = 1.5;

            for (int i = 1; i < n + 2; i++)
            {
                for (int j = 1; j <= m + 1; j++)
                {
                    dataGridView1[i, j].Value = v[i - 1, j - 1];
                }
            }

            for (int i = 1; i < n1; i++)
            {
                for (int j = 1; j < m; j++)
                {
                    if (t < Math.Abs(v[i, j] - U[i,j]))
                    {
                        t = Math.Abs(v[i, j] - U[i, j]);
                        xx = x[i];
                        yy = y[j];
                    }
                }
            }

            for (int i = n1; i < n; i++)
            {
                for (int j = 1; j < m1; j++)
                {
                    if (t < Math.Abs(v[i, j] - U[i,j]))
                    {
                        t = Math.Abs(v[i, j] - U[i, j]);
                    }
                }
            }
        }

        private void button3_Click(object sender, EventArgs e)
        {
            Form2 f = new Form2(this);
            f.Show();
        }
    }
}