import { useState, useEffect } from 'react'
import { Background } from './components/Background'
import { Card, CardHeader, CardTitle, CardContent, CardDescription, CardFooter } from './components/ui/card'
import { Button } from './components/ui/button'
import { LineChart, Line, XAxis, YAxis, CartesianGrid, Tooltip, ResponsiveContainer, ReferenceLine } from 'recharts'

interface DataPoint {
  x: number;
  u: number;
}

interface SimulationMetrics {
  duration?: string;
  stepsPerSec?: number;
  initialMass?: number;
  finalMass?: number;
  massRelError?: number;
  energyRetention?: number;
  troubledCells?: number;
  totalDgCells?: number;
  waveRegularity?: string;
}

function App() {
  const [data, setData] = useState<DataPoint[]>([])
  const [loading, setLoading] = useState(false)
  const [error, setError] = useState<string | null>(null)
  const [metrics, setMetrics] = useState<SimulationMetrics | null>(null)

  // Simulation parameter state
  const [equation, setEquation] = useState<'advection' | 'burgers'>('burgers')
  const [integrator, setIntegrator] = useState<'ssprk3' | 'euler'>('ssprk3')
  const [threads, setThreads] = useState<number>(4)

  const runSimulation = async () => {
    setLoading(true)
    setError(null)

    try {
      // Trigger simulation with configured parameters
      const runRes = await fetch('http://localhost:3001/api/run', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ equation, integrator, threads })
      })

      if (!runRes.ok) {
        const errJson = await runRes.json().catch(() => ({}))
        throw new Error(errJson.error || 'Simulation run failed')
      }

      const runJson = await runRes.json()
      if (runJson.metrics) {
        setMetrics(runJson.metrics)
      }

      // Fetch results
      const resRes = await fetch('http://localhost:3001/api/results')
      if (!resRes.ok) throw new Error('Failed to fetch simulation results')

      const jsonData = await resRes.json()
      setData(jsonData)
    } catch (err: any) {
      setError(err.message || 'An unexpected error occurred')
    } finally {
      setLoading(false)
    }
  }

  // Load existing results on startup if present
  useEffect(() => {
    fetch('http://localhost:3001/api/results')
      .then(res => res.ok ? res.json() : [])
      .then(d => {
        if (d && Array.isArray(d) && d.length > 0) setData(d)
      })
      .catch(() => {})
  }, [])

  return (
    <div className="relative min-h-screen flex items-center justify-center p-4">
      <Background />

      <main className="z-10 w-full max-w-5xl flex flex-col items-center">
        <Card className="w-full bg-slate-900/80 backdrop-blur-xl border border-white/10 text-white shadow-2xl rounded-2xl overflow-hidden">
          <CardHeader className="border-b border-white/10 pb-6">
            <div className="flex flex-col md:flex-row md:items-center md:justify-between gap-4">
              <div>
                <CardTitle className="text-3xl font-extrabold bg-clip-text text-transparent bg-gradient-to-r from-blue-400 via-indigo-400 to-purple-400">
                  HyFlow1D Web Dashboard
                </CardTitle>
                <CardDescription className="text-slate-400 mt-1">
                  High-Performance Hybrid Finite Volume / Discontinuous Galerkin Flow Testbed
                </CardDescription>
              </div>

              {/* Status Badge */}
              <div className="flex items-center gap-2">
                <span className="text-xs px-2.5 py-1 rounded-full bg-blue-500/20 text-blue-300 border border-blue-500/30 font-medium">
                  MPI + OpenMP
                </span>
                <span className="text-xs px-2.5 py-1 rounded-full bg-purple-500/20 text-purple-300 border border-purple-500/30 font-medium">
                  SSPRK3
                </span>
              </div>
            </div>
          </CardHeader>

          <CardContent className="space-y-6 pt-6">
            {/* Simulation Parameter Controls */}
            <div className="grid grid-cols-1 md:grid-cols-3 gap-4 bg-white/5 p-4 rounded-xl border border-white/10">
              <div>
                <label className="block text-xs font-semibold text-slate-300 uppercase tracking-wider mb-2">
                  Governing Equation
                </label>
                <select
                  value={equation}
                  onChange={(e) => setEquation(e.target.value as any)}
                  disabled={loading}
                  aria-label="Governing Equation"
                  className="w-full bg-slate-800 border border-white/20 rounded-lg px-3 py-2 text-sm text-white focus:outline-none focus:ring-2 focus:ring-indigo-500"
                >
                  <option value="burgers">Nonlinear Burgers' (Shock Wave)</option>
                  <option value="advection">Linear Advection (Smooth Wave)</option>
                </select>
              </div>

              <div>
                <label className="block text-xs font-semibold text-slate-300 uppercase tracking-wider mb-2">
                  Time Integrator
                </label>
                <select
                  value={integrator}
                  onChange={(e) => setIntegrator(e.target.value as any)}
                  disabled={loading}
                  aria-label="Time Integrator"
                  className="w-full bg-slate-800 border border-white/20 rounded-lg px-3 py-2 text-sm text-white focus:outline-none focus:ring-2 focus:ring-indigo-500"
                >
                  <option value="ssprk3">SSPRK3 (3rd-Order SSP RK)</option>
                  <option value="euler">Forward Euler (1st-Order)</option>
                </select>
              </div>

              <div>
                <label className="block text-xs font-semibold text-slate-300 uppercase tracking-wider mb-2">
                  OpenMP Threads: {threads}
                </label>
                <input
                  type="range"
                  min="1"
                  max="16"
                  value={threads}
                  onChange={(e) => setThreads(parseInt(e.target.value, 10))}
                  disabled={loading}
                  aria-label="OpenMP Threads"
                  className="w-full accent-indigo-500 mt-2 cursor-pointer"
                />
              </div>
            </div>

            {error && (
              <div className="bg-red-500/10 border border-red-500/30 text-red-400 px-4 py-3 rounded-lg text-sm text-center font-medium">
                Error: {error}
              </div>
            )}

            {/* Diagnostic Metrics Cards */}
            {metrics && (
              <div className="grid grid-cols-2 md:grid-cols-4 gap-3">
                <div className="bg-white/5 border border-white/10 rounded-xl p-3">
                  <div className="text-xs text-slate-400">Compute Time</div>
                  <div className="text-lg font-bold text-white mt-0.5">{metrics.duration || 'N/A'}</div>
                  <div className="text-[11px] text-slate-400">{metrics.stepsPerSec ? `${metrics.stepsPerSec.toLocaleString()} steps/s` : ''}</div>
                </div>

                <div className="bg-white/5 border border-white/10 rounded-xl p-3">
                  <div className="text-xs text-slate-400">Mass Rel Error</div>
                  <div className="text-lg font-bold text-emerald-400 mt-0.5">
                    {metrics.massRelError !== undefined ? `${metrics.massRelError.toFixed(4)}%` : 'N/A'}
                  </div>
                  <div className="text-[11px] text-slate-400">Conservative LLF</div>
                </div>

                <div className="bg-white/5 border border-white/10 rounded-xl p-3">
                  <div className="text-xs text-slate-400">Energy Retention</div>
                  <div className="text-lg font-bold text-cyan-400 mt-0.5">
                    {metrics.energyRetention !== undefined ? `${metrics.energyRetention.toFixed(2)}%` : 'N/A'}
                  </div>
                  <div className="text-[11px] text-slate-400">L² Kinetic Norm</div>
                </div>

                <div className="bg-white/5 border border-white/10 rounded-xl p-3">
                  <div className="text-xs text-slate-400">Modal Shock Sensor</div>
                  <div className={`text-lg font-bold mt-0.5 ${metrics.troubledCells && metrics.troubledCells > 0 ? 'text-amber-400' : 'text-emerald-400'}`}>
                    {metrics.troubledCells !== undefined ? `${metrics.troubledCells} / ${metrics.totalDgCells || 10}` : 'N/A'}
                  </div>
                  <div className="text-[11px] text-slate-400 truncate">
                    {metrics.waveRegularity || 'Persson-Peraire'}
                  </div>
                </div>
              </div>
            )}

            {/* Plot Canvas */}
            <div className="w-full h-[420px] bg-slate-950/60 rounded-xl border border-white/10 p-4 relative">
              <div className="absolute top-6 right-6 z-10 flex items-center gap-4 text-xs font-medium text-slate-400 bg-slate-900/80 px-3 py-1.5 rounded-lg border border-white/10 backdrop-blur-md">
                <div className="flex items-center gap-1.5">
                  <div className="w-2.5 h-2.5 rounded-full bg-blue-500/80"></div>
                  <span>Finite Volume [0, 0.5]</span>
                </div>
                <div className="flex items-center gap-1.5">
                  <div className="w-2.5 h-2.5 rounded-full bg-indigo-500"></div>
                  <span>Discontinuous Galerkin (P=3) [0.5, 1.0]</span>
                </div>
              </div>

              {data.length > 0 ? (
                <ResponsiveContainer width="100%" height="100%">
                  <LineChart data={data} margin={{ top: 20, right: 20, bottom: 20, left: 10 }}>
                    <CartesianGrid strokeDasharray="3 3" stroke="#ffffff15" />
                    <XAxis
                      dataKey="x"
                      stroke="#94a3b8"
                      tickFormatter={(val) => val.toFixed(2)}
                      domain={[0, 1]}
                    />
                    <YAxis stroke="#94a3b8" domain={['auto', 'auto']} />
                    <Tooltip
                      contentStyle={{
                        backgroundColor: 'rgba(15, 23, 42, 0.95)',
                        borderColor: 'rgba(255, 255, 255, 0.15)',
                        borderRadius: '0.5rem',
                        boxShadow: '0 10px 15px -3px rgba(0, 0, 0, 0.5)',
                      }}
                      labelStyle={{ color: '#e2e8f0', fontWeight: 'bold' }}
                      itemStyle={{ color: '#818cf8' }}
                      formatter={(val: any) => [typeof val === 'number' ? val.toFixed(5) : val, 'Solution u']}
                      labelFormatter={(label) => `Coordinate x: ${typeof label === 'number' ? label.toFixed(4) : label}`}
                    />
                    {/* Interface dividing FV and DG domains */}
                    <ReferenceLine x={0.5} stroke="#38bdf8" strokeDasharray="4 4" label={{ value: 'Interface x=0.5', fill: '#38bdf8', position: 'top' }} />
                    <Line
                      type="monotone"
                      dataKey="u"
                      stroke="#818cf8"
                      strokeWidth={2.5}
                      dot={false}
                      activeDot={{ r: 5, fill: '#c084fc' }}
                    />
                  </LineChart>
                </ResponsiveContainer>
              ) : (
                <div className="flex flex-col items-center justify-center h-full text-slate-400">
                  <p className="text-base">{loading ? 'Executing solver...' : 'No simulation data available.'}</p>
                  <p className="text-xs text-slate-500 mt-1">Configure parameters and click "Run Simulation" below.</p>
                </div>
              )}
            </div>
          </CardContent>

          <CardFooter className="flex justify-center border-t border-white/10 pt-4 pb-6">
            <Button
              onClick={runSimulation}
              disabled={loading}
              className="w-full md:w-1/3 bg-gradient-to-r from-indigo-500 to-purple-600 hover:from-indigo-600 hover:to-purple-700 text-white font-semibold py-2.5 rounded-xl shadow-lg shadow-indigo-500/25 transition-all duration-200"
            >
              {loading ? 'Simulating...' : 'Run Simulation'}
            </Button>
          </CardFooter>
        </Card>
      </main>
    </div>
  )
}

export default App
