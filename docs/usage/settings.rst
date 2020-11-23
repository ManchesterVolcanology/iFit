Settings
########

Here is an example configuration file ``example_config.yml``:

.. code-block::

  a_k: '0.0'
  a_k_fit: true
  a_w: '0.0'
  a_w_fit: true
  bgpoly_params:
  - - 0.0
    - true
  - - 0.0
    - true
  - - 0.0
    - true
  - - 1.0
    - true
  coadds: 10
  dark_flag: true
  dark_fnames: Example/dark.txt
  despike_flag: false
  fit_hi: 320.0
  fit_lo: 310.0
  flat_flag: false
  flat_path: ''
  frs_path: Ref/sao2010.txt
  fwem: '0.6'
  fwem_fit: true
  gas_params:
  - - SO2
    - 1.0e+16
    - true
    - Ref/SO2_295K.txt
  - - O3
    - 1.0e+19
    - true
    - Ref/O3_243K.txt
  - - Ring
    - 0.1
    - true
    - Ref/Ring.txt
  graph_flag: true
  graph_param: SO2
  hi_int_limit: 70000
  ils_mode: Manual
  ils_path: ''
  int_time: 100
  interp_method: linear
  k: '2.0'
  k_fit: true
  lo_int_limit: 0
  model_padding: 1.0
  model_spacing: 0.01
  ndarks: 1
  offset_params:
  - - 0.0
    - true
  resid_limit: 10.0
  resid_type: Percentage
  rt_save_path: ''
  save_path: Results/example_output.csv
  scroll_amt: 100
  scroll_flag: true
  shift_params:
  - - 0.0
    - true
  - - 0.1
    - true
  spec_fnames: 'Example/spectrum_00320.txt'
  spec_type: iFit
  spike_limit: 1000
  stray_flag: true
  stray_hi: 290.0
  stray_lo: 280.0
  update_flag: true
  wl_calib: ''
