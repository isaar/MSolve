using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Logging.Interfaces;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.Logging
{
    public class ForcesLog_v2 : IAnalyzerLog_v2
    {
        private readonly Element_v2[] elements;
        private readonly Dictionary<int, double[]> forces = new Dictionary<int, double[]>();

        public ForcesLog_v2(Element_v2[] elements)
        {
            this.elements = elements;
        }

        public Dictionary<int, double[]> Forces { get { return forces; } }
        public DateTime StartTime { get; set; }
        public DateTime EndTime { get; set; }

        public override string ToString()
        {
            StringBuilder s = new StringBuilder();
            foreach (int id in forces.Keys)
            {
                s.Append(String.Format("({0}): ", id));
                for (int i = 0; i < forces[id].Length; i++)
                    s.Append(String.Format("{0:0.00000}/", forces[id][i]));
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
                forces[e.ID] = e.ElementType.CalculateForcesForLogging(e, localVector);

                //for (int i = 0; i < stresses[e.ID].Length; i++)
                //    Debug.Write(stresses[e.ID][i]);
                //Debug.WriteLine("");
            }
        }

        #endregion
    }
}
