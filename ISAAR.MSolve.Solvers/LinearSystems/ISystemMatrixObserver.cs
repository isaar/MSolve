namespace ISAAR.MSolve.Solvers.LinearSystems
{
    /// <summary>
    /// Objects implementing this interface will be notifying before <see cref="ILinearSystem.Matrix"/> is modified.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public interface ISystemMatrixObserver
    {
        /// <summary>
        /// It will be called before setting <see cref="ILinearSystem.Matrix"/>.
        /// </summary>
        void HandleMatrixWillBeSet();
    }
}
