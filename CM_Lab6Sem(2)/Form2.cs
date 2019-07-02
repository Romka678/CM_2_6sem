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
    public partial class Form2 : Form
    {
        Form1 form;
        public Form2(Form1 form1)
        {
            InitializeComponent();
            form = form1;
            if (form.o == 0)
            {
                label3.Text = "На решение затрачено N = " + form.s1;
                label1.Text = Convert.ToString(form.eps);
                label2.Text = Convert.ToString(form.eps_max1);
                label4.Text = "Тестовая задача решена с точностью : ";
                label8.Text = Convert.ToString(form.t);
                label5.Text = "";
                pictureBox5.Visible = false;
                pictureBox6.Visible = false;
                label6.Text = "";
                label9.Text = "В точке [ " + form.xx + " ; " + form.yy + " ]";
                form.s1 = 0; form.s2 = 0; form.eps_max1 = 0; form.eps_max2 = 0;
            }
            else
            {
                label3.Text = "На решение затрачено N = " + form.s1;
                label1.Text = Convert.ToString(form.eps);
                label2.Text = Convert.ToString(form.eps_max1);
                label4.Text = "Максимальная разность двух приближений : ";
                label8.Text = Convert.ToString(form.t);
                label5.Text = label5.Text + form.s2;
                label6.Text = Convert.ToString(form.eps_max2);
                label9.Text = "В точке [ " + form.xx + " ; " + form.yy + " ]";
                form.s1 = 0; form.s2 = 0; form.eps_max1 = 0; form.eps_max2 = 0;
            }
        }
    }
}
