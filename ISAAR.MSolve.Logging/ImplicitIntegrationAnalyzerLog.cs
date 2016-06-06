using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Logging.Interfaces;

namespace ISAAR.MSolve.Logging
{
    public class ImplicitIntegrationAnalyzerLog
    {
        private List<IAnalyzerLog> logs = new List<IAnalyzerLog>();
        private readonly HashSet<Type> logTypes = new HashSet<Type>();

        public List<IAnalyzerLog> Logs { get { return logs; } }
        public void StoreResults(DateTime startTime, DateTime endTime, IAnalyzerLog log)
        {
            if (logTypes.Contains(log.GetType()) == false)
                logTypes.Add(log.GetType());
            logs.Add(log);
        }

        public void ClearResults()
        {
            logs.Clear();
        }

        public void GroupLogsByKind()
        {
            int interleaving = logTypes.Count;
            int items = logs.Count / interleaving;
            var deInterlacedLogs = new IAnalyzerLog[logs.Count];
            for (int i = 0; i < logs.Count; i++)
            {
                int item = i / interleaving;
                int group = i % interleaving;
                deInterlacedLogs[group * items + item] = logs[i];
            }
            logs = deInterlacedLogs.ToList<IAnalyzerLog>();
        }
    }
}
