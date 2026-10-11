# HyFlow1D Web Dashboard

Interactive web interface and visualization dashboard for the HyFlow1D hybrid Finite Volume / Discontinuous Galerkin aerodynamic flow solver.

## Overview

The web dashboard allows configuring and executing simulations directly through a modern web UI, rendering interactive real-time solution curves, physical conservation diagnostics, and a Three.js 3D background.

## Technology Stack

- **Runtime**: [Node.js](https://nodejs.org/) (>= 20) or [Bun](https://bun.sh/)
- **Frontend**: [React 19](https://react.dev/), [TypeScript](https://www.typescriptlang.org/)
- **Bundler**: [Vite 8](https://vitejs.dev/)
- **Styling**: [Tailwind CSS](https://tailwindcss.com/)
- **Plotting**: [Recharts](https://recharts.org/)
- **3D Graphics**: [Three.js](https://threejs.org/) via [@react-three/fiber](https://docs.pmnd.rs/react-three-fiber) and [@react-three/drei](https://github.com/pmndrs/drei)

## Architecture

- **`server.ts`**: Universal Node.js / Bun API server listening on port 3001.
  - `POST /api/run`: Spawns the compiled C++ solver (`../build/hyflow1d`) with configurable equation (`advection` or `burgers`), integrator (`ssprk3` or `euler`), and OpenMP thread count. Parses execution diagnostics and returns structured JSON metrics.
  - `GET /api/results`: Reads the generated CSV output (`solution_final.csv`) and serves coordinate points to the frontend.
- **`src/App.tsx`**: Main dashboard interface featuring parameter controls, real-time diagnostic cards (mass error, kinetic retention, Persson-Peraire modal shock sensor), and Recharts line visualizations with domain partition markers.
- **`src/components/Background.tsx`**: Three.js particle wave canvas background.

## Getting Started

### Prerequisites

Ensure the C++ solver is built:
```bash
mkdir -p build && cd build
cmake .. && make -j$(nproc)
cd ../webapp
```

### Installation

Install dependencies using `pnpm`, `npm`, or `bun`:
```bash
pnpm install
# or: npm install
# or: bun install
```

### Running the Backend Server

Start the API server (supports Node.js directly or Bun):
```bash
node server.ts
# or: bun run server.ts
```
The API server listens on `http://localhost:3001`.

### Running the Frontend

In a separate terminal, launch the Vite development server:
```bash
pnpm run dev
# or: npm run dev
# or: bun run dev
```
Open your browser at the displayed local URL (typically `http://localhost:5173`).

### Production Build

```bash
pnpm run build
pnpm run preview
```
