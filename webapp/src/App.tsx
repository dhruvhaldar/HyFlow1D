import { useState, useEffect } from 'react'
import { Background } from './components/Background'
import { Card, CardHeader, CardTitle, CardContent, CardDescription, CardFooter } from './components/ui/card'
import { Button } from './components/ui/button'
import { LineChart, Line, XAxis, YAxis, CartesianGrid, Tooltip, ResponsiveContainer } from 'recharts'

interface DataPoint {
  x: number;
  u: number;
}

function App() {
  const [data, setData] = useState<DataPoint[]>([])
  const [loading, setLoading] = useState(false)
  const [error, setError] = useState<string | null>(null)

  const runSimulation = async () => {
    setLoading(true)
    setError(null)
    setData([])

    try {
      // Trigger simulation
      const runRes = await fetch('http://localhost:3001/api/run', { method: 'POST' })
      if (!runRes.ok) throw new Error('Simulation run failed')

      // Fetch results
      const resRes = await fetch('http://localhost:3001/api/results')
      if (!resRes.ok) throw new Error('Failed to fetch results')

      const jsonData = await resRes.json()
      setData(jsonData)
    } catch (err: any) {
      setError(err.message || 'An error occurred')
    } finally {
      setLoading(false)
    }
  }

  // Optional: Load initial data if it exists
  useEffect(() => {
    fetch('http://localhost:3001/api/results')
      .then(res => res.ok ? res.json() : [])
      .then(d => {
         if (d && d.length > 0) setData(d)
      })
      .catch(() => {})
  }, [])

  return (
    <div className="relative min-h-screen flex items-center justify-center p-4">
      <Background />

      <main className="z-10 w-full max-w-4xl flex flex-col items-center">
        <Card className="w-full bg-white/10 backdrop-blur-lg border-white/20 text-white shadow-2xl">
          <CardHeader>
            <CardTitle className="text-3xl font-bold bg-clip-text text-transparent bg-gradient-to-r from-purple-400 to-pink-600">
              HyFlow1D Web Dashboard
            </CardTitle>
            <CardDescription className="text-gray-300">
              Hybrid FV-DG Simulation with 3D Background
            </CardDescription>
          </CardHeader>
          <CardContent className="space-y-4">
            <p className="text-sm leading-relaxed text-gray-200 text-center">
              Execute the C++ 1D aerodynamic flow solver directly from the web UI.
              Results are plotted in real-time using Recharts.
            </p>

            {error && <div className="text-red-400 text-center font-bold">Error: {error}</div>}

            <div className="w-full h-[400px] bg-white/5 rounded-xl border border-white/10 p-4">
              {data.length > 0 ? (
                <ResponsiveContainer width="100%" height="100%">
                  <LineChart data={data} margin={{ top: 5, right: 20, bottom: 5, left: 0 }}>
                    <CartesianGrid strokeDasharray="3 3" stroke="#ffffff33" />
                    <XAxis dataKey="x" stroke="#fff" tickFormatter={(val) => val.toFixed(2)} />
                    <YAxis stroke="#fff" />
                    <Tooltip
                      contentStyle={{ backgroundColor: 'rgba(0,0,0,0.8)', border: '1px solid rgba(255,255,255,0.2)' }}
                      labelStyle={{ color: '#fff' }}
                      itemStyle={{ color: '#8352FD' }}
                    />
                    <Line type="monotone" dataKey="u" stroke="#8352FD" strokeWidth={3} dot={false} />
                  </LineChart>
                </ResponsiveContainer>
              ) : (
                <div className="flex items-center justify-center h-full text-gray-400">
                  {loading ? 'Running simulation...' : 'No data. Run the simulation to view results.'}
                </div>
              )}
            </div>

          </CardContent>
          <CardFooter className="flex justify-center">
            <Button
              onClick={runSimulation}
              disabled={loading}
              className="w-1/2 bg-white/20 hover:bg-white/30 text-white border border-white/30 backdrop-blur-md transition-all duration-300 ease-in-out"
            >
              {loading ? 'Running...' : 'Run Simulation'}
            </Button>
          </CardFooter>
        </Card>
      </main>
    </div>
  )
}

export default App
