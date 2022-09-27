namespace optimGA 
{
    public partial class Form1 : Form
    {
        CancellationTokenSource cts = new CancellationTokenSource();

        public Form1()
        {
            InitializeComponent();
        }

        private double XSquared(double[] Values)
        {
            return (Values[0] * Values[0]);
        }

        private double Polynome(double[] Values)
        {
            return ((2 * Values[0] * Values[0])
                  + (Values[0] * Values[1])
                  + (1.5 * Values[1] * Values[1])
                  - Values[0]
                  - (2 * Values[1])
                   );
        }

        private double MaxProb(double[] Values)
        {
            double x = Values[0];
            double y = Values[1];

            return (Math.Pow(0.7 - x, 2)
                   * Math.Exp(-Math.Pow(x, 2) - Math.Pow(y + 1.3, 2))
                   - ((x - Math.Pow(x + 1, 2) - Math.Pow(y, 3))
                   * Math.Exp(-Math.Pow(x, 2) - Math.Pow(y - 0.5, 2)))
                   );
        }

        private async void button1_Click(object sender, EventArgs e)
        {
            Func<double[], double> LocalFuncPointer = XSquared;

            double[] Domain = new double[2] { -10, 10 };

            optimGA MyOptimGA = new optimGA();
            Progress<ProgressReportModel> progress = new Progress<ProgressReportModel>();
            progress.ProgressChanged += ReportProgress;
            try
            {
                double[] RetVals = await MyOptimGA.OptimGA(LocalFuncPointer, 1, Domain, progress, cts.Token);

                string Result = "";
                foreach (double r in RetVals)
                {
                    Result = Result + r.ToString() + ", ";
                }
                Result = Result.Substring(0, Result.Length - 2);
                MessageBox.Show("Expected: 0, found: " + Result + ".", "Result", MessageBoxButtons.OK, MessageBoxIcon.Information);
            }
            catch (OperationCanceledException)
            {
                System.Diagnostics.Debug.WriteLine("Operation was cancelled.");
            }
        }

        private async void button2_Click(object sender, EventArgs e)
        {
            Func<double[], double> LocalFuncPointer = Polynome;

            double[] Domain = new double[2] { -10, 10 };

            optimGA MyOptimGA = new optimGA();
            Progress<ProgressReportModel> progress = new Progress<ProgressReportModel>();
            progress.ProgressChanged += ReportProgress;
            try
            {
                double[] RetVals = await MyOptimGA.OptimGA(LocalFuncPointer, 2, Domain, progress, cts.Token);

                string Result = "";
                foreach (double r in RetVals)
                {
                    Result = Result + r.ToString() + ", ";
                }
                Result = Result.Substring(0, Result.Length - 2);
                MessageBox.Show("Expected: 0,09090909, 0,63636364, found: " + Result + ".", "Result", MessageBoxButtons.OK, MessageBoxIcon.Information);

            }
            catch (OperationCanceledException)
            {
                System.Diagnostics.Debug.WriteLine("Operation was cancelled.");
            }
        }

        private async void button3_Click(object sender, EventArgs e)
        {
            Func<double[], double> LocalFuncPointer = MaxProb;

            double[] Domain = new double[2] { -3, 3 };

            optimGA MyOptimGA = new optimGA(1);
            Progress<ProgressReportModel> progress = new Progress<ProgressReportModel>();
            progress.ProgressChanged += ReportProgress;
            try
            {
                double[] RetVals = await MyOptimGA.OptimGA(LocalFuncPointer, 2, Domain, progress, cts.Token, Eps: 0, Minimize: false);

                string Result = "";
                foreach (double r in RetVals)
                {
                    Result = Result + r.ToString() + ", ";
                }
                Result = Result.Substring(0, Result.Length - 2);
                MessageBox.Show("Expected: 0,2378043, 1,2168496, found: " + Result + ".", "Result", MessageBoxButtons.OK, MessageBoxIcon.Information);
            }
            catch (OperationCanceledException)
            {
                System.Diagnostics.Debug.WriteLine("Operation was cancelled.");
            }
        }

        private void button5_Click(object sender, EventArgs e)
        {
            cts.Cancel();
        }

        private void ReportProgress(object sender, ProgressReportModel e)
        {
            System.Diagnostics.Debug.WriteLine(Math.Round(e.Percentage, 1).ToString() + " % elapsed; " + e.SecondsElapsed + " seconds elapsed; " + e.GenerationsElapsed + " generations elapsed.");
        }
    }
}