import { Background } from './components/Background'
import { Card, CardHeader, CardTitle, CardContent, CardDescription, CardFooter } from './components/ui/card'
import { Button } from './components/ui/button'

function App() {
  return (
    <div className="relative min-h-screen flex items-center justify-center p-4">
      <Background />

      <main className="z-10 w-full max-w-md">
        <Card className="bg-white/10 backdrop-blur-lg border-white/20 text-white shadow-2xl">
          <CardHeader>
            <CardTitle className="text-3xl font-bold bg-clip-text text-transparent bg-gradient-to-r from-purple-400 to-pink-600">
              Glassmorphism UI
            </CardTitle>
            <CardDescription className="text-gray-300">
              Interactive 3D Background with React Three Fiber
            </CardDescription>
          </CardHeader>
          <CardContent className="space-y-4">
            <p className="text-sm leading-relaxed text-gray-200">
              This app features a visually stunning 3D background combined with a
              sleek glassmorphism overlay, built using Vite, Bun, TailwindCSS, and shadcn/ui.
            </p>
            <div className="flex flex-col space-y-2">
              <div className="flex items-center justify-between p-3 rounded-lg bg-white/5 border border-white/10">
                <span>Performance</span>
                <span className="font-mono text-purple-400">High</span>
              </div>
              <div className="flex items-center justify-between p-3 rounded-lg bg-white/5 border border-white/10">
                <span>Deployability</span>
                <span className="font-mono text-pink-400">Easy</span>
              </div>
            </div>
          </CardContent>
          <CardFooter>
            <Button className="w-full bg-white/20 hover:bg-white/30 text-white border border-white/30 backdrop-blur-md transition-all duration-300 ease-in-out">
              Get Started
            </Button>
          </CardFooter>
        </Card>
      </main>
    </div>
  )
}

export default App
