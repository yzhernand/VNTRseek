requires 'perl', '>= 5.24.0';
requires 'Try::Tiny';
requires 'IPC::Run';
requires 'DBD::SQLite', '>= 1.58';
requires 'Parallel::ForkManager', '>= 2.01';
on 'test' => sub {
  requires 'Test::More', '>= 1.001014';
};