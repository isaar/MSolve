using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Logging.Interfaces;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.Logging
{
    public class StressesLog_v2 : IAnalyzerLog_v2
    {
        private readonly Element_v2[] elements;
        private readonly Dictionary<int, double[]> strains = new Dictionary<int, double[]>();
        private readonly Dictionary<int, double[]> stresses = new Dictionary<int, double[]>();

        public StressesLog_v2(Element_v2[] elements)
        {
            this.elements = elements;
        }

        public Dictionary<int, double[]> Strains { get { return strains; } }
        public Dictionary<int, double[]> Stresses { get { return stresses; } }
        public DateTime StartTime { get; set; }
        public DateTime EndTime { get; set; }

        public override string ToString()
        {
            StringBuilder s = new StringBuilder();
            foreach (int id in stresses.Keys)
            {
                s.Append(String.Format("({0}): ", id));
                for (int i = 0; i < strains[id].Length; i++)
                    s.Append(String.Format("{0:e}/", strains[id][i]));
                for (int i = 0; i < stresses[id].Length; i++)
                    s.Append(String.Format("{0:0.00000}/", stresses[id][i]));
                s.Append("; ");
            }
            return s.ToString();
        }

        #region IResultStorage Members

        public void StoreResults(DateTime startTime, DateTime endTime, IVectorView solutionVector)
        {
            StartTime = startTime;
            EndTime = endTime;
            //double[] solution = ((Vector<double>)solutionVector).Data;
            foreach (Element_v2 e in elements)
            {
                double[] localVector = e.Subdomain.FreeDofOrdering.ExtractVectorElementFromSubdomain(e, solutionVector);
                var strainStresses = e.ElementType.CalculateStresses(e, localVector, 
                    new double[e.ElementType.GetElementDOFTypes(e).SelectMany(x => x).Count()]);
                strains[e.ID] = new double[strainStresses.Item1.Length];
                stresses[e.ID] = new double[strainStresses.Item2.Length];
                Array.Copy(strainStresses.Item1, strains[e.ID], strains[e.ID].Length);
                Array.Copy(strainStresses.Item2, stresses[e.ID], stresses[e.ID].Length);

                //for (int i = 0; i < stresses[e.ID].Length; i++)
                //    Debug.Write(stresses[e.ID][i]);
                //Debug.WriteLine("");
            }
        }

        #endregion
    }
}
