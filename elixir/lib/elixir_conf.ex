defmodule ElixirConf do
  @moduledoc """
  Documentation for `ElixirConf`.
  """

  @doc """
  Hello world.

  ## Examples

      iex> ElixirConf.hello()
      :world

  """

import Matrex

def hello do
    :world
end

def inner_product(v1, v2) do
  Enum.reduce(Enum.zip_with([v1, v2], fn [x1,x2] -> x1 * x2 end), &+/2)
 #Matrex.dot(v1,v2) |> (&Matrex.at(&1,1,1)).()
end

def matrix_multiply(a, b) do
  m = a[:rows]
  p = b[:cols]
  Matrex.new(m, p, fn i,j -> inner_product(row(a,i),column(b,j)) end)
end

def matrix_multiply_T(a, b) do
  bT = transpose(b)
  m = a[:rows]
  p = bT[:rows]
  Matrex.new(m, p, fn i,j -> inner_product(row(a,i),column(bT,j)) end)
end

def matrix_multiply_parallel_1(a, b) do
  m = a[:rows]
  p = b[:cols]
  parent = self()
  for i <- 1..m, j <- 1..p, do:
      spawn(fn -> send(parent, {i, j, inner_product(row(a,i),column(b,j))}) end)
  Matrex.new(m, p, fn i,j -> receive do {x,y,c} when x==i and y==j -> c end end)
end

def matrix_multiply_parallel_2(a, b) do
  m = a[:rows]
  p = b[:cols]
  parent = self()
  for i <- 1..m, do: spawn(fn -> send(parent, {i, for j <- 1..p do inner_product(row(a,i),column(b,j)) end }) end)
  Matrex.new(for i <- 1..m do receive do {x, row} when x==i -> row end end)
end


def matrix_multiply_parallel_3(a, b) do
  procs = 4
  m = a[:rows]
  p = b[:cols]
  parent = self()
  for i0 <- Enum.map(0..procs-1, fn x -> x*div(m,procs)+1 end), do:
      spawn(fn -> send(parent, {i0, calculate_stripe_C(a, b, p, i0, div(m,procs))}) end)
  Matrex.new(collect_C(div(m,procs), procs))
end

def calculate_stripe_C(a, b, p, i0, s) do
  for i <- i0..i0+s-1 do
      for j <- 1..p do
          inner_product(row(a,i), column(b,j))
      end
  end
end

def collect_C(ms, p) do
   List.foldl(for k <- 0..p-1 do
                  receive do {x, stripe} when x == k*ms+1 -> stripe end
              end,
              [], fn (x,acc) -> x ++ acc end)
end

#def recv_stripe_3(_,  _, 0), do: []
#def recv_stripe_3(i, ms, p), do: (receive do {ii, stripe} when ii==i-> stripe end) ++ recv_stripe_3(i + ms, ms, p - 1)

def matrix_multiply_parallel_4(a, b) do
  procs_x = 2
  procs_y = 2
  m = a[:rows]
  p = b[:cols]
  parent = self()
  for i0 <- Enum.map(0..procs_x-1, fn x -> div(x*m,procs_x)+1 end), j0 <- Enum.map(0..procs_y-1, fn x -> div(x*p,procs_y)+1 end), do:
      spawn(fn -> for i <- i0..i0+div(m,procs_x)-1, j <- j0..j0+div(p,procs_y)-1, do: send(parent, {i, j,  inner_product(row(a,i),column(b,j))}) end)
  Matrex.new(m, p, fn i,j -> receive do {ii,jj,cij} when ii==i and jj==j -> cij end end)
end

# Cannon's algorithm

# WORKER

def shift(x0, pid) do
    send({:worker, pid}, x0)
    receive do x1 -> x1 end
end

def cannon_skewing(x0, _, 0), do: x0
def cannon_skewing(x0, pid, i) do
    x1 = shift(x0, pid)
    cannon_skewing(x1, pid, i-1)
end

def cannon_shifting(_, _, c0, _, _, 0), do: c0
def cannon_shifting(a0, b0, c0, pidU, pidL, i) do
    b1 = shift(b0, pidU)
    a1 = shift(a0, pidL)
    c1 = Matrex.add(c0, matrix_multiply(a1,b1))
    cannon_shifting(a1, b1, c1, pidU, pidL, i-1)
end

def cannon(i, j, a0, b0, c0, pidU, pidL, p) do
    b1 = cannon_skewing(b0, pidU, i-1)
    a1 = cannon_skewing(a0, pidL, j-1)
    cannon_shifting(a1, b1, c0, pidU, pidL, p)
end

def matrix_multiply_cannon_worker(i, j, pidM, pidU, pidL, m, p) do
  Process.register(self(), :worker)
  a = Matrex.ones(m)
  b = Matrex.ones(m)
  c0 = Matrex.zeros(m)
  c1 = cannon(i, j, a, b, c0, pidU, pidL, p)
  x = c1[1][1]
  send({:manager, pidM}, x)
end


#MANAGER

def matrix_multiply_cannon_manager(nodes, m, u, l) do
  Process.register(self(), :manager)
  procs = length(nodes)
  p = round(:math.sqrt(procs))
  create_mesh_row(nodes, u, l, m, p, p)
  receive_block(procs)
end

def receive_block(0), do: []
def receive_block(p) do
  c = receive do c0 -> c0 end
  [c | receive_block(p-1)]
end

def create_mesh_row(_, _, _, _, 0, _), do: 0
def create_mesh_row(nodes0, u, l, m, i, p) do
  nodes1 = create_mesh_col(nodes0, u, l, m, p, i, p)
  create_mesh_row(nodes1, u, l, m, i-1, p)
end

def create_mesh_col(node, _, _, _, _, _, 0), do: node
def create_mesh_col([node|nodes], u, l, m, p, i, j) do
  Node.spawn_link(node, ElixirConf, :matrix_multiply_cannon_worker, [i, j, node(), u.(node), l.(node), m, p])
  create_mesh_col(nodes, u, l, m, p, i, j-1)
end

# WORKPOOL

def matrix_multiply_parallel_5(a, b) do
  procs = 4
  m = a[:rows]
  p = b[:cols]
  pIdC = spawn(ElixirConf, :collect_C, [self(), m, p])

  for _ <- 1..procs do
    spawn(ElixirConf, :matrix_multiply_workpool_worker, [pIdC, self()])
  end

  for i <- 1..m do
    for j <- 1..p do
      as = row(a,i)
      bs = column(b,j)
      receive do {:request, pidW} -> send(pidW, {:job, as, bs, i, j}) end
    end
  end

  for _ <- 1..procs do
    receive do {:request, pidW} -> send(pidW, :finish) end
  end

  receive do {:final, c} -> c end
end

def matrix_multiply_workpool_worker(pidC, pidM) do
  send(pidM, {:request, self()})
  receive do
      {:job, a, b, i, j} -> send(pidC, {:result, i, j, inner_product(a,b)})
                            matrix_multiply_workpool_worker(pidC, pidM)
      :finish -> :ok
  end
end

def collect_C(pidM, m, p) do
  x = Matrex.new(m, p, fn i,j -> receive do {:result, x, y, c} when x==i and y==j -> c end end)
  send(pidM, {:final, x})
end

end
