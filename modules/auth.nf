nextflow.enable.dsl=2

process GISAID_AUTHENTICATION {
    output:
    path "auth_complete.txt", emit: auth_status

    script:
    """
    python << 'PYTHON_SCRIPT'
    import sys
    from outbreak_data import authenticate_user
    try:
        authenticate_user.get_authentication()
        print('Authentication already exists')
    except:
        print('Missing authentication. Please run `python gisaid_authentication.py` to authenticate and try again')
        sys.exit(1)
    PYTHON_SCRIPT
    
    if [ \$? -eq 0 ]; then
        echo "Authentication verified" > auth_complete.txt
    else
        exit 1
    fi
    """
}

workflow AUTH {
    main:
    auth_status = GISAID_AUTHENTICATION()

    emit:
    auth_status
}