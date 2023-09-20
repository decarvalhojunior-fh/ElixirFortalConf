%%%-------------------------------------------------------------------
%% @doc erlangconf public API
%% @end
%%%-------------------------------------------------------------------

-module(erlangconf_app).

-behaviour(application).

-export([start/2, stop/1, matrix_multiply/2, matrix_multiply_T/2, matrix_multiply_parallel_1/2, 
         matrix_multiply_parallel_2/2, matrix_multiply_parallel_3/2, 
         matrix_multiply_cannon_manager/4, matrix_multiply_cannon_worker/7, cannon_skewing_2/4, shift_2/3, matrix_multiply_parallel_5/2,
         matrix_multiply_workpool_worker/2, collect_C_2/3, collect_C/3, adaptative_quadrature/4, adaptative_quadrature_worker_request/3]).

start(_StartType, _StartArgs) ->
    erlangconf_sup:start_link().

stop(_State) -> ok.

%% internal functions

inner_product(V1, V2) -> numerl:vec_dot(V1,V2).

matrix_multiply_T(A,BT) -> 
    {_,M,_,_} = A,
    {_,P,_,_} = BT,
    numerl:matrix([[inner_product(numerl:row(I,A), numerl:row(J,BT)) || J <- lists:seq(1,P) ] || I <- lists:seq(1,M)]).

matrix_multiply(A,B) -> 
    {_,M,_,_} = A,
    {_,_,P,_} = B,
    numerl:matrix([[inner_product(numerl:row(I,A), numerl:col(J,B)) || J <- lists:seq(1,P) ] || I <- lists:seq(1,M)]).


matrix_multiply_parallel_1(A, B) -> 
    {_,M,_,_} = A,
    {_,__,P,_} = B,
    Pid = self(),
    create_tasks_I_1(A, B, M, P, Pid),
    numerl:matrix(fetch_result_I_1(M,P)).

create_tasks_I_1(_, _, I, _, _) when I==0 -> ok;
create_tasks_I_1(A, B, I, P, Pid) when I>0 -> 
    create_tasks_J_1(A, B, I, P, Pid),
    create_tasks_I_1(A, B, I-1, P, Pid).

create_tasks_J_1(_, _, _, J, _) when J==0 -> ok;
create_tasks_J_1(A, B, I, J, Pid) when J>0 ->
    spawn(fun () ->  Pid ! {I, J, inner_product(numerl:row(I,A), numerl:col(J,B))} end),
    create_tasks_J_1(A, B, I, J-1, Pid).

fetch_result_I_1(I, _) when I==0 -> [];
fetch_result_I_1(I, P) when I>0 -> 
    Row = fetch_result_J_1(I, P), [Row | fetch_result_I_1(I-1, P)].

fetch_result_J_1(_, J) when J==0 -> [];
fetch_result_J_1(I, J) when J>0 -> 
    receive {X,Y,C} when X==I andalso Y==J -> [C | fetch_result_J_1(I, J-1)] end.

fetch_results_1(M, P) -> [[receive {X,Y,C} when X==I andalso Y==J -> C end || I <- lists:seq(1,M)] || J <- lists:seq(1,P)].


matrix_multiply_parallel_2(A, B) -> 
    {_,M,_,_} = A,
    {_,_,P,_} = B,
    Pid = self(),
    create_tasks_I_2(A, B, M, P, Pid),
    numerl:matrix([receive {X,C} when X==I -> C end || I <- lists:seq(1,M)]).

create_tasks_I_2(_, _, I, _, _) when I==0 -> ok;
create_tasks_I_2(A, B, I, P, Pid) when I>0 -> 
    spawn(fun () -> V = [inner_product(numerl:row(I,A), numerl:col(J,B)) || J <- lists:seq(1,P)], Pid ! {I, V} end),
    create_tasks_I_2(A, B, I-1, P, Pid).

matrix_multiply_parallel_3(A, B) -> 
    PROCS = 4,
    {_,M,_,_} = A,
    {_,_,P,_} = B,
    Pid = self(),
    create_tasks_S(A, B, 1, M div PROCS, P, Pid, PROCS),
    numerl:matrix(collect_C(M div PROCS, PROCS)).

create_tasks_S(_, _, _, _, _, _, PROCS) when PROCS==0 -> ok;
create_tasks_S(A, B, I0, Ms, P, Pid, PROCS) when PROCS>0 -> 
    spawn(fun () -> S = calculate_stripe_C(A, B, I0, Ms, P), Pid ! {I0, S} end),
    create_tasks_S(A, B, I0+Ms, Ms, P, Pid, PROCS-1).

calculate_stripe_C(A, B, I0, Ms, P) -> 
    [[inner_product(numerl:row(I,A), numerl:col(J,B)) 
           || J <- lists:seq(1,P)] 
                     || I <- lists:seq(I0,I0+Ms-1)].

collect_C(Ms, PROCS) -> lists:foldl(fun (X,ACC) -> X ++ ACC end, [], [receive {X, C} when X == K*Ms+1 -> C end || K <- lists:seq(0,PROCS-1)]).


% Canon's Algorithm

% WORKERS (mesh processes)

shift(X0, Pid) -> 
    {worker, Pid} ! X0, 
    receive  X1 -> X1 end.

shift_2(X0, Pid1, Pid2) -> Pid2 ! shift(X0, Pid1).

cannon_skewing(X0, _, I) when I==0 -> X0;
cannon_skewing(X0, Pid, I) when I>0 ->
    X1 = shift(X0, Pid),
    cannon_skewing(X1, Pid, I-1).

cannon_skewing_2(X0, Pid1, I, Pid2) -> Pid2 ! cannon_skewing(X0, Pid1, I).

cannon_shifting(_, _, C0, _, _, I) when I==0 -> C0;
cannon_shifting(A0, B0, C0, PidU, PidL, I) when I > 0-> 
 %   spawn(erlangconf_app, shift_2, [B0, PidU, self()]),
 %   spawn(erlangconf_app, shift_2, [A0, PidL, self()]),
 %   B1 = receive B -> B end,
 %   A1 = receive A -> A end,
    B1 = shift(B0, PidU),
    A1 = shift(A0, PidL),
    C1 = numerl:add(C0, matrix_multiply(A1,B1)),
    cannon_shifting(A1, B1, C1, PidU, PidL, I-1).

cannon(I, J, A0, B0, C0, PidU, PidL, P) ->
 %   spawn(erlangconf_app, cannon_skewing_2, [B0, PidU, P, self()]),
 %   spawn(erlangconf_app, cannon_skewing_2, [A0, PidL, P, self()]),
 %   B1 = receive B -> B end,
 %   A1 = receive A -> A end,
    B1 = cannon_skewing(B0, PidU, I-1),
    A1 = cannon_skewing(A0, PidL, J-1),
    cannon_shifting(A1, B1, C0, PidU, PidL, P).

matrix_multiply_cannon_worker (I, J, PidM, PidU, PidL, M, P)  -> 
    register(worker, self()),
    A = numerl:matrix([[1 || Y <- lists:seq(1,M)] || Y <- lists:seq(1,M)]), %numerl:rnd_matrix(M),
    B = numerl:matrix([[1 || Y <- lists:seq(1,M)] || Y <- lists:seq(1,M)]), %numerl:rnd_matrix(M),
    C0 = numerl:zeros(M,M), 
    C1 = cannon(I, J, A, B, C0, PidU, PidL, P),
    X = numerl:get(1,1,C1),
    {manager, PidM} ! X.

% MANAGER (launches the processes accross cluster nodes)

matrix_multiply_cannon_manager(Nodes, M, U, L) -> 
    register(manager, self()),
    PROCS = length(Nodes),
    P = round(math:sqrt(PROCS)),
    create_mesh_row(Nodes, U, L, M, P, P),
    receive_block(PROCS).

receive_block(P) when P==0 -> [];
receive_block(P) when P>0 -> [receive C0 -> C0 end | receive_block(P-1)].

create_mesh_row(_, _, _, _, I, _) when I==0 -> ok;
create_mesh_row(Nodes0, U, L, M, I, P) when I>0 ->
    Nodes1 = create_mesh_col(Nodes0, U, L, M, P, I, P),
    create_mesh_row(Nodes1, U, L, M, I-1, P).

create_mesh_col(Node, _, _, _, _, _, J) when J==0 -> Node;
create_mesh_col([Node|Nodes], U, L, M, P, I, J) when J>0 ->
    spawn(Node, erlangconf_app, matrix_multiply_cannon_worker, [I, J, node(), U(Node), L(Node), M, P]),
    create_mesh_col(Nodes, U, L, M, P, I, J-1).

% WORKPOOL


matrix_multiply_parallel_5(A, B) ->
    PROCS = 4, 
    {_,M,_,_} = A,
    {_,_,P,_} = B,
    PIdC = spawn(erlangconf_app, collect_C, [self(), M, P]),
    launch_workers(PIdC, PROCS),
    workpool_conversation_loop_i(A, B, 1, 1, M, P, PROCS),
    receive {mm, C} -> C end.    

launch_workers(_, K) when K==0 -> ok;
launch_workers(PIdC, K) when K>0 ->
    spawn(erlangconf_app, matrix_multiply_workpool_worker, [PIdC, self()]),
    launch_workers(PIdC, K-1).

workpool_conversation_loop_i(_, _, I, _, M, _, PROCS) when I==M+1 -> finish_workers(PROCS);
workpool_conversation_loop_i(A, B, I, J, M, P, PROCS) when I=<M ->
    workpool_conversation_loop_j(A, B, I, 1, P),
    workpool_conversation_loop_i(A, B, I+1, J, M, P, PROCS).

workpool_conversation_loop_j(_, _, _, J, P) when J==P+1 -> ok;
workpool_conversation_loop_j(A, B, I, J, P) when J=<P -> 
    As = numerl:row(I,A),
    Bs = numerl:col(J,B),
    receive
        {request, PidW} -> PidW ! {job, As, Bs, I, J}, 
                           workpool_conversation_loop_j(A, B, I, J+1, P)
    end.

finish_workers(K) when K==0 -> ok;
finish_workers(K) when K>0 -> 
    receive
        {request, PidW} -> PidW ! finish,
                           finish_workers(K-1) 
    end.

matrix_multiply_workpool_worker(PidC, PidM) ->
    PidM ! {request, self()},
    receive 
        {job, A, B, I, J} -> PidC ! {result, I, J, inner_product(A,B)}, 
                             matrix_multiply_workpool_worker(PidC, PidM);
        finish -> ok  
    end.

collect_C(PidM, M, P) -> PidM ! {mm, numerl:matrix(collect_C_I(1,M,P))}.

collect_C_I(I, M,  _) when I==M+1 -> [];
collect_C_I(I, M,  P) when I=<M -> 
    Row = collect_C_IJ(I, 1, P), [Row | collect_C_I(I+1, M, P)].

collect_C_IJ(_, J, P) when J==P+1 -> [];
collect_C_IJ(I, J, P) when J=<P ->
     receive {result,X,Y,C} when X==I andalso Y==J -> [C | collect_C_IJ(I, J+1, P)] end.

collect_C_2(PidM, M, P) -> PidM ! {mm, numerl:matrix([[receive {result,X,Y,C} when X==I andalso Y==J -> C end || I <- lists:seq(1,M)] || J <- lists:seq(1,P)])}.


adaptative_quadrature(F, L, R, EPSILON) ->
      PROCS = 4,
      launch_aqworkers(PROCS, F, EPSILON),
      FL = F(L), 
      FR = F(R),
      A = ((FL - FR) * (R-L)) / 2,
      workpool_conversation_loop_aq([{job, L, R, FL, FR, A}], 0.0).

workpool_conversation_loop_aq([], TOTAL) -> TOTAL;
workpool_conversation_loop_aq([J|JS], TOTAL) ->
     receive 
        {request, WPid} -> WPid ! J, workpool_conversation_loop_aq(JS, TOTAL);
        {job, L, R, FL, FR, A} -> workpool_conversation_loop_aq([{job, L, R, FL, FR, A} | [J | JS]], TOTAL);
        {total, S} -> workpool_conversation_loop_aq(JS, TOTAL + S)
     end.
      
launch_aqworkers(P, _, _) when P==0 -> ok;
launch_aqworkers(P, F, EPSILON) when P>0 -> 
    spawn(erlangconf_app, adaptative_quadrature_worker_request, [self(), F, EPSILON]),
    launch_aqworkers(P-1, F, EPSILON).

adaptative_quadrature_worker_request(Parent, F, EPSILON) -> 
         Parent ! {request, self()},
         receive {job, L, R, FL, FR, A} ->
              adaptative_quadrature_worker_continue(Parent, F, L, R, FL, FR, A, EPSILON)
         end.

adaptative_quadrature_worker_continue(Parent, F, L, R, FL, FR, LRA, EPSILON) -> 
              M =  L+R div 2,
              FM = F(M),
              LA = ((FL - FM) * (M-L)) / 2,
              RA = ((FM - FR) * (R-M)) / 2,
              if 
                abs((LA + RA) - LRA) > EPSILON -> 
                    parent ! {job, L, M, FL, FM, LA},
                    adaptative_quadrature_worker_continue(Parent, F, M, R, FM, FR, RA, EPSILON);
                true -> Parent ! {total, LRA},
                    adaptative_quadrature_worker_request(Parent, F, EPSILON)
              end.
