using System;
using System.Xml.Linq;

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

        private double rndRange(double Min, double Max)
        {
            Random rnd = new Random();
            return ((rnd.NextDouble() * (Max - Min)) + Min);
        }

        private double rndGauss(double Mu, double Sigma, double Bounds)
        {
            double x = 0;
            double y = 0;
            double Density = 0;
            Random rnd = new Random();

            do
            {
                x = (2 * Bounds * rnd.NextDouble()) + (Mu - Bounds);
                y = rnd.NextDouble();
                Density = (1 / Math.Sqrt(2 * Math.PI * Sigma)) * Math.Exp(-(Math.Pow(x - Mu, 2)) / (2 * Sigma));
            } while (y < Density);

            return (x);
        }

        private void QuickSort(ref double[,] vSort, int Index = 0, int Start = -1, int End = -1)
        {
            if (Start == -1) Start = 0;
            if (End == -1) End = vSort.GetLength(0) - 1;

            int i;
            int j;
            double h;
            double x;

            i = Start;
            j = End;
            x = vSort[(Start + End) / 2, Index];

            do
            {
                while (vSort[i, Index] < x) i++;
                while (vSort[j, Index] > x) j--;

                if (i <= j)
                {
                    for (int u = 0; u < vSort.GetLength(1); u++)
                    {
                        h = vSort[i, u];
                        vSort[i, u] = vSort[j, u];
                        vSort[j, u] = h;
                    }
                    i++;
                    j--;
                }
            } while (i <= j);

            if (Start < j) QuickSort(ref vSort, Index, Start, j);
            if (i < End) QuickSort(ref vSort, Index, i, End);
        }

        private double[,] ReverseArray(double[,] vRev)
        {
            double[,] Target = new double[vRev.GetLength(0), vRev.GetLength(1)];

            for (int i = 0; i < vRev.GetLength(0); i++)
            {
                for (int j = 0; j < vRev.GetLength(1); j++)
                {
                    Target[i, j] = vRev[vRev.GetLength(0) - 1 - i, j];
                }
            }

            return (Target);
        }

        private int[] Resample(int Num)
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

        private async Task<double[]> Evaluate(double[,] Population, Func<double[], double> Fun)
        {
            int PopSize = Population.GetLength(0);
            int nVars = Population.GetLength(1) - 1;

            List<Task<double>> tasks = new List<Task<double>>();

            for (int i = 0; i < PopSize; i++)
            {
                double[] Row = new double[nVars];
                for (int j = 0; j < nVars; j++)
                {
                    Row[j] = Population[i, j];
                }
                tasks.Add(Task.Run(() => Fun(Row)));
            }

            double[] RetVals = await Task.WhenAll(tasks);

            return(RetVals);
        }

        // Public methods ------------------------------------------------------------------

        public async Task<double[]> OptimGA(Func<double[], double> Fun, int nVars, double[] Domain, IProgress<ProgressReportModel> progress, CancellationToken cancellationToken, int PopSize = 100, int MaxGenerations = 1000, double MutationProbability = 0.7, double RecombineProbability = 1.0, double Eps = 0.00001, bool Minimize = true)
        {
            bool Max = !Minimize;
            double[,] Population = new double[PopSize, nVars + 1];
            double[,] TempPopulation = new double[PopSize, nVars + 1];
            double[] CurrVars = new double[nVars];
            double[] retVal = new double[nVars];
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
                for (int j = 0; j < nVars; j++)
                {
                    Population[i, j] = rndRange(Domain[0], Domain[1]);
                }
            }

            for (int i = 0; i < MaxGenerations; i++)
            {
                if (_timeBudget != 0)
                {
                    if (DateTime.Now.Subtract(StartTime).Seconds > TimeBudget)
                    {
                        for (int j = 0; j < nVars; j++)
                        {
                            retVal[j] = Population[0, j];
                        }

                        _endTime = DateTime.Now;
                        _iters = i;
                        _terminationReason = TerminationReason.TimeOut;
                        return (retVal);
                    }
                }

                double[] EvRes = await Evaluate(Population, Fun);
                for (int j = 0; j < PopSize; j++) Population[j, nVars] = EvRes[j];

                QuickSort(ref Population, nVars);
                if (Max) Population = ReverseArray(Population);

                report.Percentage = i / MaxGenerations * 100;
                report.SecondsElapsed = DateTime.Now.Subtract(StartTime).Seconds;
                report.GenerationsElapsed = i;
                progress.Report(report);

                cancellationToken.ThrowIfCancellationRequested();

                EpsEmp = Math.Abs(Population[0, nVars] - Population[(PopSize / 2) - 1, nVars]);
                if (EpsEmp < Eps)
                {
                    for (int j = 0; j < nVars; j++)
                    {
                        retVal[j] = Population[0, j];
                    }
                    _endTime = DateTime.Now;
                    _iters = i;
                    _terminationReason = TerminationReason.ConvergenceReached;
                    return (retVal);
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
                            TempPopulation[j + (PopSize / 2), k] = ((Population[RecombineSet1[j], k] + Population[RecombineSet2[j], k]) / 2) + (Mutate * rndGauss(0, Domain[1] - Domain[0], 5));
                            if (TempPopulation[j + (PopSize / 2), k] < Domain[0]) TempPopulation[j + (PopSize / 2), k] = Domain[0];
                            if (TempPopulation[j + (PopSize / 2), k] > Domain[1]) TempPopulation[j + (PopSize / 2), k] = Domain[1];
                        }
                    }
                    else
                    {
                        for (int k = 0; k < nVars; k++)
                        {
                            TempPopulation[j + (PopSize / 2), k] = Population[j, k] + (Mutate * rndGauss(0, Domain[1] - Domain[0], 5));
                            if (TempPopulation[j + (PopSize / 2), k] < Domain[0]) TempPopulation[j + (PopSize / 2), k] = Domain[0];
                            if (TempPopulation[j + (PopSize / 2), k] > Domain[1]) TempPopulation[j + (PopSize / 2), k] = Domain[1];
                        }
                    }
                }
                Population = TempPopulation;
                _iters = i;
            }

            for (int j = 0; j < nVars; j++)
            {
                retVal[j] = Population[0, j];
            }
            _endTime = DateTime.Now;
            _terminationReason = TerminationReason.MaxItersReached;
            return (retVal);
        }
    }
}
