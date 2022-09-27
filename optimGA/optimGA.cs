namespace optimGA
{
    public enum TerminationReason
    {
        ConvergenceReached,
        TimeOut,
        MaxItersReached
    }

    public class ProgressReportModel
    {
        public double Percentage { get; set; }
        public int SecondsElapsed { get; set; }
        public int GenerationsElapsed { get; set; }
    }

    public class optimGA
    {
        // Internal variables --------------------------------------------------------------

        private int _timeBudget = 300; // Measured in seconds.
        private DateTime _startTime;
        private DateTime _endTime;
        private int _iters;
        private TerminationReason _terminationReason;

        // Getters and setters -------------------------------------------------------------

        public int TimeBudget
        {
            get
            {
                return (_timeBudget);
            }
            set
            {
                _timeBudget = value;
            }
        }

        public DateTime StartTime
        {
            get
            {
                return (_startTime);
            }
        }

        public DateTime EndTime
        {
            get
            {
                return (_endTime);
            }
        }

        public int ElapsedSeconds
        {
            get
            {
                return (_endTime.Subtract(_startTime).Seconds);
            }
        }

        public int NumberOfIterations
        {
            get
            {
                return (_iters);
            }
        }

        public TerminationReason ReasonOfTermination
        {
            get
            {
                return (_terminationReason);
            }
        }

        // Constructors --------------------------------------------------------------------

        // Parameterless constructor.
        public optimGA()
        {
        
        }

        public optimGA(int TimeBudgetInSeconds)
        {
            _timeBudget = TimeBudgetInSeconds;
        }

        // Private methods -----------------------------------------------------------------

        private static double rndRange(double Min, double Max)
        {
            Random rnd = new Random();
            return ((rnd.NextDouble() * (Max - Min)) + Min);
        }

        private static double rndGauss(double Mu, double Sigma, double Bounds)
        {
            double x;
            double y;
            double Density;
            Random rnd = new Random();

            do
            {
                x = (2 * Bounds * rnd.NextDouble()) + (Mu - Bounds);
                y = rnd.NextDouble();
                Density = (1 / Math.Sqrt(2 * Math.PI * Sigma)) * Math.Exp(-(Math.Pow(x - Mu, 2)) / (2 * Sigma));
            } while (y < Density);

            return (x);
        }

        private static int[] Resample(int Num)
        {
            int[] OutVec = new int[Num];
            Random rnd = new Random();
            int Candidate = rnd.Next(0, Num);

            Array.Fill(OutVec, -1);

            for (int i = 0; i < Num; i++)
            {
                while (Array.Exists(OutVec, element => element == Candidate))
                {
                    Candidate = rnd.Next(0, Num);
                }
                OutVec[i] = Candidate;
            }

            return (OutVec);
        }

        private static async Task<double[]> Evaluate(List<double[]> Population, Func<double[], double> Fun)
        {
            int PopSize = Population.Count;
            int nVars = Population[0].Length - 1;

            List<Task<double>> tasks = new List<Task<double>>();

            foreach (double[] r in Population)
            {
                double[] Row = r;
                Array.Resize(ref Row, Row.Length - 1);
                tasks.Add(Task.Run(() => Fun(Row)));
            }

            double[] RetVals = await Task.WhenAll(tasks);

            return(RetVals);
        }

        // Public methods ------------------------------------------------------------------

        public async Task<double[]> OptimGA(Func<double[], double> Fun, int nVars, double[] Domain, IProgress<ProgressReportModel> progress, CancellationToken cancellationToken, int PopSize = 100, int MaxGenerations = 1000, double MutationProbability = 0.7, double RecombineProbability = 1.0, double Eps = 0.00001, bool Minimize = true)
        {
            List<double[]> Population = new List<double[]>();
            List<double[]> TempPopulation = new List<double[]>();
            double[] Row;
            double EpsEmp;
            int[] RecombineSet1;
            int[] RecombineSet2;
            double Mutate;
            double Recombine;
            Random rnd = new Random();

            ProgressReportModel report = new ProgressReportModel();

            _startTime = DateTime.Now;

            for (int i = 0; i < PopSize; i++)
            {
                Row = new double[nVars + 1];
                for (int j = 0; j < nVars; j++)
                {
                    Row[j] = rndRange(Domain[0], Domain[1]);
                }
                Population.Add(Row);
            }

            for (int i = 0; i < MaxGenerations; i++)
            {
                if (_timeBudget != 0)
                {
                    if (DateTime.Now.Subtract(StartTime).Seconds > TimeBudget)
                    {
                        Row = Population[0];
                        Array.Resize(ref Row, Row.Length - 1);

                        _endTime = DateTime.Now;
                        _iters = i;
                        _terminationReason = TerminationReason.TimeOut;
                        return (Row);
                    }
                }

                double[] EvRes = await Evaluate(Population, Fun);
                for (int j = 0; j < PopSize; j++) Population[j][nVars] = EvRes[j];

                Population.Sort((double[] a, double[] b) => a.Last().CompareTo(b.Last()));
                if (!Minimize) Population.Reverse();

                report.Percentage = ((double) i / MaxGenerations) * 100;
                report.SecondsElapsed = DateTime.Now.Subtract(StartTime).Seconds;
                report.GenerationsElapsed = i;
                progress.Report(report);

                cancellationToken.ThrowIfCancellationRequested();

                EpsEmp = Math.Abs(Population[0][nVars] - Population[(PopSize / 2) - 1][nVars]);
                if (EpsEmp < Eps)
                {
                    Row = Population[0];
                    Array.Resize(ref Row, Row.Length - 1);

                    _endTime = DateTime.Now;
                    _iters = i;
                    _terminationReason = TerminationReason.ConvergenceReached;
                    return (Row);
                }

                RecombineSet1 = Resample(PopSize / 2);
                RecombineSet2 = Resample(PopSize / 2);

                TempPopulation = Population;
                for (int j = 0; j < PopSize / 2; j++)
                {
                    Mutate = rnd.NextDouble();
                    if (Mutate >= MutationProbability) Mutate = 1; else Mutate = 0;

                    Recombine = rnd.NextDouble();
                    if (Recombine >= (1 - RecombineProbability))
                    {
                        for (int k = 0; k < nVars; k++)
                        {
                            TempPopulation[j + (PopSize / 2)][k] = ((Population[RecombineSet1[j]][k] + Population[RecombineSet2[j]][k]) / 2) + (Mutate * rndGauss(0, Domain[1] - Domain[0], 5));
                            if (TempPopulation[j + (PopSize / 2)][k] < Domain[0]) TempPopulation[j + (PopSize / 2)][k] = Domain[0];
                            if (TempPopulation[j + (PopSize / 2)][k] > Domain[1]) TempPopulation[j + (PopSize / 2)][k] = Domain[1];
                        }
                    }
                    else
                    {
                        for (int k = 0; k < nVars; k++)
                        {
                            TempPopulation[j + (PopSize / 2)][k] = Population[j][k] + (Mutate * rndGauss(0, Domain[1] - Domain[0], 5));
                            if (TempPopulation[j + (PopSize / 2)][k] < Domain[0]) TempPopulation[j + (PopSize / 2)][k] = Domain[0];
                            if (TempPopulation[j + (PopSize / 2)][k] > Domain[1]) TempPopulation[j + (PopSize / 2)][k] = Domain[1];
                        }
                    }
                }
                Population = TempPopulation;
                _iters = i;
            }

            Row = Population[0];
            Array.Resize(ref Row, Row.Length - 1);

            _endTime = DateTime.Now;
            _terminationReason = TerminationReason.MaxItersReached;
            return (Row);
        }
    }
}
