# How to contribute
1. Copy branch name from linear. It will look like `zpp-2137-do-something`.
2. Create branch locally with `git checkout -b zpp-2137-do-something`. 
3. Code. 
4. Create branch with given name on origin with `git push --set-upstream origin zpp-2137-do-something`. 
5. Create a pull request to master. 


# Environment variables
MMseq uses env varibales to connect to db: 
`"PG_HOST",
"PG_PORT",
"PG_DBNAME",
"PG_USER",
"PG_PASSWORD"`
