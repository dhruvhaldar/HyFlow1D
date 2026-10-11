# HyFlow1D Web Dashboard

Interactive web interface and visualization dashboard for the HyFlow1D hybrid Finite Volume / Discontinuous Galerkin aerodynamic flow solver.

## Overview

The web dashboard allows running simulations directly through a modern web UI, rendering interactive real-time solution curves and a Three.js 3D background.

## Technology Stack

- **Runtime**: [Bun](https://bun.sh/)
- **Frontend**: [React 19](https://react.dev/), [TypeScript](https://www.typescriptlang.org/)
- **Bundler**: [Vite 8](https://vitejs.dev/)
- **Styling**: [Tailwind CSS](https://tailwindcss.com/)
- **Plotting**: [Recharts](https://recharts.org/)
- **3D Graphics**: [Three.js](https://threejs.org/) via [@react-three/fiber](https://docs.pmnd.rs/react-three-fiber) and [@react-three/drei](https://github.com/pmndrs/drei)

## Architecture

- **`server.ts`**: Lightweight Bun API server running on port 3001.
  - `POST /api/run`: Triggers the compiled C++ solver (`../build/hyflow1d --clean`).
  - `GET /api/results`: Reads the generated CSV output (`solution_final.csv`) and serves structured JSON data points to the frontend.
- **`src/App.tsx`**: Main dashboard interface featuring execution controls, real-time status indicators, and Recharts line visualizations.
- **`src/components/Background.tsx`**: Three.js particle wave canvas background.

## Getting Started

### Prerequisites

Ensure [Bun](https://bun.sh/) is installed:
```bash
curl -fsSL https://bun.sh/install | bash
```

Ensure the C++ solver is built:
```bash
cd ..
mkdir -p build && cd build
cmake .. && make
cd ../webapp
```

### Installation

Install frontend and server dependencies:
```bash
bun install
```

### Running the Backend Server

Start the Bun API server:
```bash
bun run server.ts
```
The API server listens on `http://localhost:3001`.

### Running the Frontend

In a separate terminal, launch the Vite development server:
```bash
bun run dev
```
Open your browser at the displayed local URL (typically `http://localhost:5173`).

### Production Build

```bash
bun run build
bun run preview
```
